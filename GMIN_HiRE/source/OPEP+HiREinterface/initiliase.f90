!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Adapted for use with GMIN and OPTIM from original
! kr366> reading parameters and initialising the potential 
!        (formerly read_parameters and md_initialise)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initialise(CONF,RNAT,DNAT,PROT)
    USE COMMONS, ONLY: MYUNIT,SAXST
    USE DEFINITIONS
    USE RNAparams, only: molecule_type
    USE GEOMETRIC_CORRECTIONS, ONLY: CENTER_OF_MASS
    USE PORFUNCS
    USE ION_PAIR, ONLY: ION_PAIR_SCALING, ION_PAIR_CONTROL
    use prec_hire  
    IMPLICIT NONE
    INTEGER :: I1, PARAMS_UNIT, GETUNIT, FCHAIN
    LOGICAL :: EXISTT
    LOGICAL, INTENT(IN) :: RNAT,DNAT,PROT

    type(t_conformations), intent(inout) :: conf 

  ! We now read the number of fragments 
    FCHAIN = GETUNIT()
    INQUIRE(FILE='ichain.dat',EXIST=EXISTT)
    IF (EXISTT) THEN
       OPEN(UNIT=FCHAIN,FILE='ichain.dat',STATUS='OLD',ACTION='READ')
    ELSE
       INQUIRE(FILE='ichain_RNA.dat',EXIST=EXISTT)
       IF (EXISTT) THEN
          OPEN(UNIT=FCHAIN,FILE='ichain_RNA.dat',STATUS='OLD',ACTION='READ')
       ELSE
          WRITE(MYUNIT,'(A)') 'initialise> ichain.dat does not exist'
          STOP
       ENDIF
    ENDIF

    READ(FCHAIN,*) NFRAG
    WRITE(MYUNIT,'(A,I6)') 'initialise> read ichain for restriction', NFRAG     
    ALLOCATE(LIST_FRAGMENTS(NFRAG,2))   
    DO I1=1, NFRAG
       READ(FCHAIN,*) LIST_FRAGMENTS(I1,1),LIST_FRAGMENTS(I1,2)
    ENDDO
    CLOSE(FCHAIN)
    VECSIZE = 3 * NATOMS
    VECSIZE1 =  VECSIZE / NFRAG

!set to defaults for now
    force_scaling_factor = 1.0D0
    PBC = .FALSE.
    BL = 1.0D0
    C_M = .FALSE.
    ion_pair_control = .FALSE.
    ion_pair_scaling = 1.0D0
    N_REPLICA = 1

    WRITE(MYUNIT,'(A39,I12)')   'initialise> Number of atoms             : ', NATOMS
    WRITE(MYUNIT,'(A39,I12)')   'initialise> Number of fragments         : ', NFRAG
    WRITE(MYUNIT,'(A39,F12.6)') 'initialise> Potential scaling factor    : ', force_scaling_factor
    WRITE(MYUNIT,'(A39,L12)')   'initialise> Periodic Boundary Condition : ', PBC
    WRITE(MYUNIT,'(A39,F12.6)') 'initialise> Box Length                  : ', BL
    WRITE(MYUNIT,'(A39,L12)')   'initialise> center of mass for pdb      : ', C_M


    !CALL read_parameters_md()
    ALLOCATE(pos(VECSIZE))       
    ALLOCATE(posref(VECSIZE))       
    ALLOCATE(force(VECSIZE))       
    ALLOCATE(atomic_type(NATOMS))
    ALLOCATE(mass(vecsize))

  ! We first set-up pointers for the x, y, z components in the position and
  ! forces

    x    => pos(1:3*natoms:3)
    y    => pos(2:3*natoms:3)
    z    => pos(3:3*natoms:3)

    xref => posref(1:3*natoms:3)
    yref => posref(2:3*natoms:3)
    zref => posref(3:3*natoms:3)

    fx   => force(1:3*natoms:3)
    fy   => force(2:3*natoms:3)
    fz   => force(3:3*natoms:3)

! We then initialise the potential
    RNA_simulation = (rnat.or.dnat)
    prot_simulation = prot
    if (rnat) molecule_type = "RNA"
    if (dnat) molecule_type = "DNA"

    call initialise_potential()
    call center_of_mass(natoms,pos,mass)
    ! We get the information on the bonds and their ideal length
    call readtop1(nbondh,nbonda,nbonds)
    nbondt = nbondh + nbonda

    ! Defines the various vectors necessary
    if (.not.allocated(ia)) allocate(ia(nbondt), ib(nbondt))
    if (.not.allocated(redu1_mass))  &
    allocate(redu1_mass(nbondt), redu2_mass(nbondt), reduA(nbondt), reduB(nbondt), reduAA(nbondt), reduBB(nbondt))
    if (.not.allocated(beq)) allocate(beq(nbondt), beq2(nbondt), ibeq2(nbondt), rij1(3,nbondt), rij2(nbondt))
    call readtop2(nbondt,ia,ib,beq,beq2) !redu1_mass,redu2_mass)
    degree_freedom = dble(vecsize) - 6.0d0

   RETURN
END SUBROUTINE initialise

subroutine initialise_potential()
    use commons, only: debug
    use definitions
    use prec_hire
    implicit none
    
    integer :: i
    real(kind = real64), dimension(vecsize) :: xpos    ! Positions ( (x1,y1,z1),(x2,y2,z2)...
    real(kind = real64), dimension(natoms)  :: amass   ! 1 over the masses
  
    ! We now call the protein part to get the positions and force

    if (prot_simulation) then
      call initialise_protein(PBC,BL,C_M,natoms,xpos,amass,atomic_type,force_scaling_factor, DEBUG)
    endif

    if (RNA_simulation) then
      call initialise_RNA(PBC,BL,C_M,natoms,xpos,amass,atomic_type,force_scaling_factor, DEBUG)
    endif

    ! We reorder the positions so that they are consistent with the rest of the program
    do i = 1, natoms
       x(i) = xpos(3*i-2)
       y(i) = xpos(3*i-1)
       z(i) = xpos(3*i)
    end do
    pos = xpos

    mass(1:3*natoms:3) = 1.0d0/amass
    mass(2:3*natoms:3) = 1.0d0/amass
    mass(3:3*natoms:3) = 1.0d0/amass
  
    mass = mass * 2390.0 ! From amu to kcal/mol fs^2 / A^2 
    return
end subroutine initialise_potential

SUBROUTINE initialise_SAXS()
    use commons, only: SAXST, SAXSPRINT, SAXSNPRINT, SAXSMODULT, SAXSINVSIG, &
        &SAXSSOLT, REFINET, NWATLAY, SAXSMAX, WATRAD
    use definitions, only: SAXSs, SAXSc, SAXS_print, compute_SAXS_serial, &
        &modulate_SAXS_serial, SAXS_invsig
    use SAXS_scoring  
    implicit none
    integer :: getunit

    compute_saxs_serial = SAXST
    saxs_print = SAXSPRINT
    modulate_saxs_serial = SAXSMODULT
    saxs_invsig = SAXSINVSIG
    SAXS_max = SAXSMAX
    in_solution_curve = SAXSSOLT
    refine_hydration_layer = REFINET
    DX = WATRAD
    n_shells = NWATLAY
    
    SAXSs = getunit()
    open(SAXSs, file='SAXS_score.dat', status='unknown',action='write',position='append')
    SAXSc = getunit()
    open(SAXSc, file='SAXS_curve.dat', status='unknown',action='write',position='append')
    call set_SAXS_scoring()
END SUBROUTINE initialise_SAXS
    
! H-bonds detection
SUBROUTINE initialise_HB()
    use commons, only: RNAHBT
    use definitions, only: do_hb, save_hb, hbdat, ihb
    do_hb=RNAHBT
    save_hb=.false.
    ihb=0
    hbdat = getunit()
    open(hbdat, file='hbonds.dat', status='unknown',action='write',position='append')
END SUBROUTINE initialise_HB

!clean-up routine for end of run
SUBROUTINE end_definitions()
    USE commons, only: SAXST, RNAHBT
    USE definitions, only: pos, posref, force, atomic_type, mass, SAXSs, SAXSc, HBdat
    USE RNAparams

    IF (SAXST) THEN
        CLOSE(SAXSS)
        CLOSE(SAXSC)
    END IF
    IF (RNAHBT) CLOSE(HBDAT)

    IF (ALLOCATED(pos)) DEALLOCATE(pos)
    IF (ALLOCATED(posref)) DEALLOCATE(posref)
    IF (ALLOCATED(force)) DEALLOCATE(force)
    IF (ALLOCATED(atomic_type)) DEALLOCATE(atomic_type)
    IF (ALLOCATED(mass)) DEALLOCATE(mass)

    DEALLOCATE(AMASS,IAC,LABRES, FATORTIT)
    DEALLOCATE(RK,REQ,TK,TEQ,PK,PN,IPN,PHASE,GAMC,GAMS,FMN)        
    DEALLOCATE(IB,JB,ICB,IT,JT,KT,ICT,IP,JP,KP,LP,ICP)

    RETURN
END SUBROUTINE end_definitions



