!> @file
!> Contain HIRE_INTERFACE module providing interface to HiRE library

!> Module to access the HiRE library, this is the interface that should be accessed from third-party programmers 
MODULE HIRE_INTERFACE
   IMPLICIT NONE
   !> Custom real, doube precision 64-bit kind
   INTEGER, PARAMETER  :: R64 = SELECTED_REAL_KIND(15, 307)

   CONTAINS
      !> Subroutine to initialise HiRE
      !> @brief
      !>
      !> Routines needs the topology and the scale.dat file names.\n
      !> The setup for HiRE is then called with these files, and the number of atoms returned.
      !> After calling this subroutine, the HiRE interface can be used.
      !>
      !> @param[in] TOPNAME - topology name
      !> @param[in] SCALEDATNAME - name fo the scale.dat file
      !> @param[out] NATOMS - number of CG particles in the module
      !>
      !> @see SETUP
      !> @see RESET_POT_ENE
      SUBROUTINE HIRE_INITIALISE(TOPNAME, SCALEDATNAME, NATOMS)
#ifdef __HIRE
         USE CALCFORCES, ONLY: EVEC, RESET_POT_ENE
         USE  VAR_DEFS, ONLY: NPARTICLES
#endif
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: TOPNAME !Name of topology file 
         CHARACTER(LEN=*), INTENT(IN) :: SCALEDATNAME !Name of scale.dat file     
         INTEGER, INTENT(OUT) :: NATOMS

         NATOMS = 0
#ifdef __HIRE       
         !call setup subroutine for HiRe
         CALL HIRE_SETUP(TOPNAME, SCALEDATNAME)
         !initialise all energies to zero
         CALL RESET_POT_ENE(EVEC)
        
         NATOMS = NPARTICLES
#endif
      END SUBROUTINE HIRE_INITIALISE 
  
      !> Routine to get energy and gradient
      !> @brief
      !>
      !> Subroutine to interface call to the calculation of the energy and gradient.\n
      !> The coordinates are not stored internally, so any new call is independent of previous ones.\n
      !>
      !> @attention The force is set to 1.0D10 and the gradient to 0.0D0. 
      !> If these values are returned something went wrong with the energy and gradient calculations!
      !>
      !> @param[in] NOPT - number of degress of freedom (required to get the array sizes correct, 
      !> should be 3*NATOMS as returned from HIRE_INITIALISE)
      !> @param[in] X - input coordinates
      !> @param[out] E - total HiRE energy for this set of coordinates
      !> @param[out] GRAD - gradient array (size is NOPT)
      !> @param[in] SAXST - use SAXs force?
      !>
      !> @see CALCFORCE
      SUBROUTINE HIRE_ENERGY_GRAD(NOPT, X, E, GRAD, SAXST)
#ifdef __HIRE
         USE CALCFORCES, ONLY: CALCFORCE
         USE SAXS_DEFS, ONLY: SAXSFORCET
#endif   
         INTEGER, INTENT(IN) :: NOPT                    !should be 3*NATOMS
         REAL(KIND = R64), INTENT(IN) :: X(NOPT)     !input coordinates
         REAL(KIND = R64), INTENT(OUT) :: E          !energy
         REAL(KIND = R64), INTENT(OUT) :: GRAD(NOPT) !gradient
         LOGICAL, INTENT(IN) :: SAXST

         GRAD(:) = 0.0D0
         E = 1.0D10     
#ifdef __HIRE       
         ! set whether we calculate SAXS
         SAXSFORCET = SAXST
         !call subroutine to calculate force and energy
         CALL CALCFORCE(NOPT,X,GRAD,E)
         GRAD(:)=-GRAD(:)
#endif
      END SUBROUTINE HIRE_ENERGY_GRAD

      !> Routine to obtain a numerical Hessian matrix
      !> @brief
      !>
      !> The subroutine provides a Hessian estimate for a given structure.\n
      !> The subroutine uses central finite difference with 4 terms (Richardson interpolation).
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] COORDS - input coordinates
      !> @param[in] DELTA - step-size taken for interpolation, the four points are 
      !> COORDS + 2*DELTA, COORDS + DELTA, COORDS - DELTA and COORDS - 2*DELTA.
      !> @param[out] HESSIAN - estimate for Hessian matrix (size is NOPT by NOPT)
      !>
      !> @see HIRE_ENERGY_GRAD
      SUBROUTINE HIRE_NUMHESS(NOPT,COORDS,DELTA,HESSIAN)  
         !input and output variables
         INTEGER, INTENT(IN) :: NOPT
         REAL(KIND = R64), INTENT(IN) :: COORDS(NOPT), DELTA
         REAL(KIND = R64), INTENT(OUT) :: HESSIAN(NOPT,NOPT)
         !variables used internally
         REAL(KIND = R64), DIMENSION(NOPT) ::  COORDS_PLUS, COORDS_PLUS2, &
                                               COORDS_MINUS, COORDS_MINUS2
         REAL(KIND = R64), DIMENSION(NOPT) :: GRAD_PLUS, GRAD_PLUS2, &
                                              GRAD_MINUS, GRAD_MINUS2
         INTEGER :: I
         REAL(KIND = R64) :: DUMMY_ENERGY 
       
         HESSIAN(:,:) = 0.0D0
         !use central finite difference with 4 terms (Richardson interpolation)
         DO I = 1,NOPT
            COORDS_PLUS(:) = COORDS(:)
            COORDS_PLUS(I) = COORDS(I) + DELTA
            CALL HIRE_ENERGY_GRAD(NOPT,COORDS_PLUS,DUMMY_ENERGY,GRAD_PLUS,.FALSE.)
            COORDS_PLUS2(:) = COORDS(:)
            COORDS_PLUS2(I) = COORDS(I) + 2.0D0 * DELTA
            CALL HIRE_ENERGY_GRAD(NOPT,COORDS_PLUS2,DUMMY_ENERGY,GRAD_PLUS2,.FALSE.)
            COORDS_MINUS(:) = COORDS(:)
            COORDS_MINUS(I) = COORDS(I) - DELTA
            CALL HIRE_ENERGY_GRAD(NOPT,COORDS_MINUS,DUMMY_ENERGY,GRAD_MINUS,.FALSE.)
            COORDS_MINUS2(:) = COORDS(:)
            COORDS_MINUS2(I) = COORDS(I) - 2.0D0 * DELTA
            CALL HIRE_ENERGY_GRAD(NOPT,COORDS_MINUS2,DUMMY_ENERGY,GRAD_MINUS2,.FALSE.)
            HESSIAN(I,:) = (GRAD_MINUS2(:) - 8.0D0 * GRAD_MINUS(:) &
                           + 8.0D0 *GRAD_PLUS(:) - GRAD_PLUS2(:))/(12.0D0*DELTA)
         END DO            
      END SUBROUTINE HIRE_NUMHESS

      !> Print energy decomposition
      !> @brief
      !>
      !> Routines calls current energies and prints them to a provided output unit.\n
      !> The energy printed are the values from the last call to CALCFORCE. Especially, 
      !> if Hessians are calculated care needs to be taken, that the correct values are reported here.
      !>
      !> @param[in] OUTUNIT - output unit
      !>
      !> @see PRINT_POT_ENE
      SUBROUTINE HIRE_DECOMPE(OUTUNIT)
#ifdef __HIRE
         USE CALCFORCES, ONLY: PRINT_POT_ENE, EVEC
#endif
         !output unit for writing
         INTEGER, INTENT(IN) :: OUTUNIT
#ifdef __HIRE       
         CALL PRINT_POT_ENE(OUTUNIT, EVEC)
#endif
      END SUBROUTINE HIRE_DECOMPE 

      !> Create pdb file
      !> @brief
      !>
      !> This subroutine creates a pdb file with the current molecular data loaded in HiRE and the provided coordinates.\n
      !> There is an option to produce chain identifiers in the pdb.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] X - input coordinates
      !> @param[in] PDBNAME - name for pdb to be saved
      !> @param[in] CHAINIDT - chain identifier to be used in pdb file?
      !>
      !> @see PDBFROMX
      SUBROUTINE DUMP_PDB(NOPT,X,PDBNAME,CHAINIDT)
#ifdef __HIRE
         USE PDB_OUT, ONLY: PDBFROMX
#endif
         INTEGER, INTENT(IN) :: NOPT
         REAL(KIND=R64), INTENT(IN) :: X(NOPT)
         CHARACTER(LEN=*), INTENT(IN) :: PDBNAME
         LOGICAL, INTENT(IN) :: CHAINIDT
#ifdef __HIRE
         CALL PDBFROMX(X,PDBNAME,CHAINIDT)
#endif
      END SUBROUTINE DUMP_PDB

      !> Routine to return CG particle names
      !>
      !> @param[in] NPART - number of CG particles
      !> @param[out] PNAMES - names of CG particles
      SUBROUTINE PASS_PARTICLE_NAMES(NPART,PNAMES)
#ifdef __HIRE
         USE VAR_DEFS, ONLY: IGRAPH
#endif
         INTEGER, INTENT(IN) :: NPART
         CHARACTER(4), INTENT(OUT) :: PNAMES(NPART)
         PNAMES(1:NPART) = ""
#ifdef __HIRE
         PNAMES(1:NPART) = IGRAPH(1:NPART)
#endif
      END SUBROUTINE PASS_PARTICLE_NAMES

      !> Routine to return CG particle masses
      !>
      !> @param[in] NPART - number of CG particles
      !> @param[out] MASSES - masses of CG particles      
      SUBROUTINE PASS_HIRE_MASSES(NPART,MASSES)
#ifdef __HIRE
      USE VAR_DEFS, ONLY: AMASS
#endif
         INTEGER, INTENT(IN) :: NPART
         REAL(KIND = R64), INTENT(OUT) :: MASSES(NPART)
         MASSES(1:NPART) = 0.0D0
#ifdef __HIRE
         MASSES(1:NPART) = AMASS(1:NPART)
#endif
      END SUBROUTINE PASS_HIRE_MASSES

      !> Subroutine to terminatew HiRE interface (deallocates all arrays in the library)
      SUBROUTINE TERMINATE_HIRE()
#ifdef __HIRE
         CALL HIRE_FINISH()
#endif
      END SUBROUTINE TERMINATE_HIRE

      !> Subroutine to initialise SAXS calculations
      SUBROUTINE SETUP_SAXS(SAXST_,SAXSPRINT_,SAXSMODULT_,SAXSINVSIG_,SAXSMAX_, &
                          SAXSSOLT_,REFINET_,WATRAD_,NWATLAY_)
#ifdef __HIRE                        
         USE MOD_INIT, ONLY: INITIALISE_SAXS
         USE SAXS_DEFS
         USE SAXS_CALCS, ONLY: SAXS_SETUP
#endif
         LOGICAL, INTENT(IN) :: SAXST_, SAXSPRINT_,SAXSMODULT_,SAXSSOLT_,REFINET_
         INTEGER, INTENT(IN) :: NWATLAY_
         REAL(KIND = R64), INTENT(IN) :: SAXSINVSIG_, WATRAD_, SAXSMAX_
#ifdef __HIRE
         compute_saxs_serial = SAXST_
         saxs_print = SAXSPRINT_
         modulate_saxs_serial = SAXSMODULT_
         saxs_invsig = SAXSINVSIG_
         SAXSMAX = SAXSMAX_
         SAXSSOLT = SAXSSOLT_
         REFINET = REFINET_
         WATRAD = WATRAD_
         NWATLAY = NWATLAY_
         CALL INITIALISE_SAXS()

#endif
      END SUBROUTINE SETUP_SAXS  
   
      ! Subroutine to pass updated variables for SAXS curve modulation
      SUBROUTINE MODULATE_SAXS(SAXSINVSIG_,SAXSWAVE_,SAXSOFFI_,SAXSMODI_)
#ifdef __HIRE
         USE SAXS_DEFS, ONLY: SAXSINVSIG,SAXSWAVE,SAXSOFFI,SAXSMODI
#endif
         REAL(KIND=R64), INTENT(IN) :: SAXSINVSIG_ , SAXSWAVE_, SAXSOFFI_, SAXSMODI_
#ifdef __HIRE
         SAXSINVSIG = SAXSINVSIG_
         SAXSWAVE = SAXSWAVE_
         SAXSOFFI  = SAXSOFFI_
         SAXSMODI = SAXSMODI_
#endif
      END SUBROUTINE MODULATE_SAXS

      !> Routine to get SAXS force
      !> @brief
      !>
      !> Subroutine to interface call to the calculation of the SAXS force.\n
      !>
      !> @param[in] NOPT - number of degress of freedom (required to get the array sizes correct, 
      !> should be 3*NATOMS as returned from HIRE_INITIALISE)
      !> @param[in] X - input coordinates
      !> @param[out] GRAD - SAXS force array (size is NOPT)
      SUBROUTINE HIRE_SAXS_FORCE(NOPT, X, ESAXS, GRAD, GRADT, HYDRATET)
#ifdef __HIRE
         USE MOD_SAXS, ONLY: RNA_SAXS_FORCE
         USE SAXS_DEFS, ONLY: SAXSFORCET
#endif   
         INTEGER, INTENT(IN) :: NOPT                    !should be 3*NATOMS
         REAL(KIND = R64), INTENT(IN) :: X(NOPT)     !input coordinates
         REAL(KIND = R64), INTENT(OUT) :: GRAD(NOPT) !gradient
         REAL(KIND = R64), INTENT(OUT) :: ESAXS     
         LOGICAL, INTENT(IN) :: GRADT   
         LOGICAL, INTENT(IN) :: HYDRATET   
         GRAD(:) = 0.0D0
         ESAXS = 1.0D10     
#ifdef __HIRE       
         ! set whether we calculate SAXS
         SAXSFORCET = .TRUE.
         !call subroutine to calculate force and energy
         CALL RNA_SAXS_FORCE(NOPT, X, ESAXS, GRAD, GRADT, HYDRATET)
         GRAD(:)=-GRAD(:)
#endif
      END SUBROUTINE HIRE_SAXS_FORCE

      !> Routine to initialise hydrogen-bonding information 
      SUBROUTINE HIRE_HB_INIT()
#ifdef __HIRE
         USE HB_DEFS
         USE UTILS_IO, ONLY: GETUNIT
       
         DO_HB = .TRUE.
         HBDAT = GETUNIT()
         OPEN(HBDAT, FILE = 'hbonds.dat', STATUS='UNKNOWN', ACTION='WRITE', POSITION = 'APPEND')
#endif
      END SUBROUTINE HIRE_HB_INIT

      !> Calculate and produce hydrogen-bonding information for inout coordinates 
      SUBROUTINE CALL_HBDAT(NOPT,X)
#ifdef __HIRE
         USE CALCFORCES, ONLY: CALC_HBONDS
#endif
         INTEGER, INTENT(IN) :: NOPT
         REAL(KIND = R64), INTENT(IN) :: X(NOPT)
#ifdef __HIRE      
         CALL CALC_HBONDS(NOPT,X)
#endif
      END SUBROUTINE CALL_HBDAT

      !> Routine to return number of dihedral angles
      SUBROUTINE GET_NDIHS(NTORS)
#ifdef __HIRE
         USE MOD_DIHEDRALS, ONLY: NDIHS
#endif
         INTEGER, INTENT(OUT) :: NTORS
         NTORS = 0
#ifdef __HIRE
         NTORS=NDIHS
#endif
      END SUBROUTINE GET_NDIHS

      !> Pass dihedral information from HiRE
      SUBROUTINE PASS_DIHINFO(NTORS,DIHINFO)
#ifdef __HIRE
         USE MOD_DIHEDRALS, ONLY: GET_DIHEDRALS
#endif
         INTEGER, INTENT(IN) :: NTORS
         INTEGER, INTENT(OUT) :: DIHINFO(NTORS,5)

         DIHINFO(1:NTORS,1:5) = 0
#ifdef __HIRE
         CALL GET_DIHEDRALS(DIHINFO)
#endif
      END SUBROUTINE PASS_DIHINFO

      SUBROUTINE SET_SCALING(SCALEVAL)
#ifdef __HIRE
         USE CALCFORCES, ONLY: SCALING
#endif     
         REAL(KIND = R64), INTENT(IN) :: SCALEVAL(8)

#ifdef __HIRE
         SCALING(1:8) = SCALEVAL(1:8)
#endif         
      END SUBROUTINE SET_SCALING

      SUBROUTINE SET_UNIV_SCALING(SCALEVAL)
#ifdef __HIRE
         USE CALCFORCES, ONLY: SCALING
#endif     
         REAL(KIND = R64), INTENT(IN) :: SCALEVAL
            
#ifdef __HIRE
         SCALING(1:8) = SCALEVAL
#endif         
      END SUBROUTINE SET_UNIV_SCALING    

END MODULE HIRE_INTERFACE
