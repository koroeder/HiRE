module definitions
  ! 
  ! kr366> adapted for use with GMIN and OPEP from original defs.F90 and md_main_module.f90
  !
  use prec_hire
  implicit none

  save
  ! Lists of parameters
  integer :: NFRAG 
  integer :: NATOMS 
  integer :: VECSIZE
  integer :: VECSIZE1                     ! Length of the force and position vectors of one fragment 

  character(len=5), dimension(:), allocatable :: atomic_type ! Atomic type
  real(kind = real64), dimension(:), allocatable, target :: force       ! Working forces on the atoms
  real(kind = real64), dimension(:), allocatable, target :: pos         ! Working positions of the atoms
  real(kind = real64), dimension(:), allocatable, target :: posref      ! Reference position
  real(kind = real64), dimension(:), allocatable, target :: mass        ! masses
  real(kind = real64) :: force_scaling_factor  ! Factor for rescaling the potential

  real(kind = real64), dimension(:), pointer :: x, y, z
  real(kind = real64), dimension(:), pointer :: xref, yref, zref
  real(kind = real64), dimension(:), pointer :: fx, fy, fz

  integer, dimension(:,:), allocatable :: list_fragments

  real(kind = real64) :: total_energy, ref_energy       ! Energies

  integer :: T_id
  integer :: E_scale
  logical :: init_single_file
  logical :: PBC  ! periodic boundary condition
  logical :: C_M  ! center of mass for writing non-broken chains in pdb file
  real(kind = real64) :: BL   ! box length 
  real(kind = real64) :: degree_freedom

  logical :: prot_simulation, RNA_simulation
  integer :: N_prot, N_RNA, N_replica
  logical :: use_tit

  integer :: flag_tit
  integer :: n_steps_tit


  !kr366> migrated from md_defs in main MD module
  integer :: nbondh, nbonda, nbondt, nbonds
  real(kind = real64), allocatable, dimension(:) :: beq, beq2, ibeq2, rij2
  real(kind = real64), allocatable, dimension(:,:) :: rij1
  real(kind = real64), allocatable, dimension(:) :: redu1_mass, redu2_mass, reduA, reduB, reduAA, reduBB
  integer, allocatable, dimension(:) :: ia, ib
  real(kind = real64) :: epspos, epsvel

  !lm759> hbonds detection
  logical :: do_hb, save_hb
  integer :: hbdat, ihb

  ! saxs serial parameters
  logical :: compute_SAXS_serial, SAXS_save, modulate_SAXS_serial
  real(kind = real64) :: SAXS_invsig
  integer :: saxs_serial_step, SAXSs, SAXSc
  integer :: n_rate_saxs

  ! Booleans deciding if we should compute a SAXS forces next step and save SAXS data
  logical :: calc_SAXS_force = .true.
  logical :: SAXS_print = .false.
  

  type t_conformations
    character(len=20) :: path
    character(len=20) :: logfile
    integer ::  id, counter
    real(kind = real64), dimension(:), allocatable :: pos, posref, vel
    real(kind = real64), dimension(:), allocatable :: temperatures,scales
    real(kind = real64) :: energy, temperature, energyscale,free_energy,E_energy,N_energy, score
  end type t_conformations

end module definitions

module numerical_defs
    use prec_hire
    implicit none
    real(kind = real64), parameter :: pi = 3.141592653589793d+00
    real(kind = real64), parameter :: rad2deg = 57.29577951308232088d0
    real(kind = real64), parameter :: rad2 = rad2deg*rad2deg

    integer MAXPRE, MAXNAT, MAXTTY, MAXXC, MAXPNB
    integer MAXBO, MAXTH, MAXPHI, MAXPAI
    parameter (MAXPRE = 1500)               !! maximum number of residues
    parameter (MAXNAT = MAXPRE*6)           !! maximum number of atoms
    parameter (MAXTTY = 50000)              !! maximum number of residue name types
    parameter (MAXXC = 3*MAXNAT)            !! maximum number of cart coord
    parameter (MAXPNB = 3*MAXPRE*MAXPRE)    !! max number of SC-SC interactions
    parameter (MAXBO  = MAXNAT)             !! maximum number of bonds
    parameter (MAXTH = MAXNAT*3)            !! maximum number of bond angles
    parameter (MAXPHI = MAXNAT*4)           !! maximum number of torsional angles
    parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)!! max number of nonbonded-pairs

    save
end module numerical_defs

module hydro_defs
    use prec_hire
    use numerical_defs, only: maxpnb
    implicit none

    real(kind = real64) :: rncoe(maxpnb), vamax(maxpnb), epshb_mcmc(maxpnb)
    integer :: ni(maxpnb),nj(maxpnb),ivi(maxpnb),ivj(maxpnb), nstep, nb

    save
end module hydro_defs

module score
  use prec_hire
    implicit none
    real(kind = real64), dimension(272) :: score_prot
    real(kind = real64), dimension(25)  :: score_RNA
    save
end module score

module PBC_defs
  use prec_hire
  implicit none
  real(kind = real64) box_length, inv_box_length
  logical periodicBC,CM
  save

  contains

!==============================================================================
!=======================    periodic boundary condition  ======================
!=======================        for each atom            ======================
!==============================================================================
!     Application of the nearest-image convention.
!     the atoms are place in a box going from
!     [-0.5*box_length, 0.5*box_length]
  elemental function pbc_mic(x) !allows any array size (function must be pure)
    implicit none
    real(kind = real64), intent(in) :: x
    real(kind = real64) :: pbc_mic

    pbc_mic = x - box_length * dnint( x * inv_box_length)
  end function

    ! same as above, 3 coordinates at a time
    subroutine pbc_mic3( x )
    implicit none

    real(kind = real64), dimension(3), intent(inout) :: x

    x = x - box_length * dnint(x * inv_box_length)

    end subroutine pbc_mic3

    ! apply same pbc as p on all posv
    subroutine pbc_ref(p, posv)
      use definitions
      implicit none

      real(kind = real64), dimension(3), intent(in) :: p
      real(kind = real64), dimension(:), intent(inout) :: posv
      integer, dimension(3) :: pr
      integer :: i

      pr = dnint(p * inv_box_length)
      do i = 1, natoms
        posv(i*3-2:i) = posv(i*3-2:i) - box_length * pr
      enddo

    end subroutine pbc_ref

end module PBC_defs

module PDBtext
    use numerical_defs
    implicit none
    character(7) :: text2(MAXNAT)
    character(5) :: text3(MAXNAT)
    character(7) :: text4(MAXNAT)
    integer :: numres(MAXNAT),Id_atom(MAXNAT)
    
    save
end module PDBtext

module fragments
    use numerical_defs
    implicit none
    integer :: nfrag, nfrag_prot, nfrag_rna,lenfrag(MAXPRE),ichain(MAXNAT)

    save
end module fragments


!================
! PROTEIN MODULES
!================
 
module system_defs_prot
  use prec_hire
    use numerical_defs
    implicit none
    ! misc1
    integer NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,&
                NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA

    real(kind = real64) amass(maxnat)
    integer iac(maxnat), nno(maxtty)

    integer IPRES(MAXPRE)
    character(4) :: IGRAPH(MAXNAT),LABRES(MAXNAT)


    !  COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC
    ! --- SET SOME PARAMETERS
    real(kind = real64) :: &
      CUT    = 100.0,&    !! NO CUTOFF
      SCNB   = 80,&       !! divide the 1-4 VDW interactions by 8.0
      SCEE   = 80.0,&     !! divide the 1-4 ELEC interactions by 2.0
      IDIEL  = 0.0,&      !! 2.0 dielectric constant epsilon = 2r
      DIELC  = 1.0        !! ...............................

    save
end module system_defs_prot

module param_defs_prot
  use prec_hire
    use numerical_defs
    implicit none

    real(kind = real64) :: &
      RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),&
      PK(MAXPHI),PN(MAXPHI),&
      PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),&
      GAMC(MAXPHI),GAMS(MAXPHI),FMN(MAXPHI)
    integer :: IPN(MAXPHI)

    integer :: &
      IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),&
      ICBH(MAXBO)

    integer :: &
      IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH),&
      JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)

    integer :: &
      IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),ICP(MAXPHI)

    integer :: &
      IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)

    integer :: NUMEX(MAXNAT),NATEX(MAXTTY)

    integer :: NUMBND,NUMANG,NPTRA,NPHB

    integer :: NPAIR2,IPAIR(MAXPAI),JPAIR(MAXPAI)

    save
end module param_defs_prot

module cutoffs_prot
  use prec_hire
    implicit none
    real(kind = real64) :: &
      rcut2_caca_scsc_out, rcut2_caca_scsc_in,&
      rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,&
      rcut2_4b_out, rcut2_4b_in,&
      rcut2_lj_out, rcut2_lj_in

    save
end module cutoffs_prot

module charge_prot
  use prec_hire
    use numerical_defs
    implicit none

    real(kind = real64) :: CG(MAXNAT)

    save
end module charge_prot


module protein_common
  use prec_hire
  use numerical_defs
  implicit none

  real(kind = real64), parameter :: zero = 0.0d0
  real(kind = real64), parameter :: one = 1.0d0
  real(kind = real64), parameter :: six = 6.0d0
  real(kind = real64), parameter :: ten = 1.0d1
  real(kind = real64), parameter :: twelve = 12.0d0
  real(kind = real64), parameter :: pt999 = 0.9990d0  
  integer :: ialpha(MAXNAT),ibeta(MAXNAT),icoeff(MAXPAI)

  real(kind = real64) :: walpha(20), wbeta(20), walpha_foal(20), wbeta_fobe(20)

  real(kind = real64) :: resp(MAXPRE), EHHB1, EHHB3

  real(kind = real64) :: ct0lj(maxpnb), ct2lj(maxpnb)

end module protein_common


