module SAXS_scoring
  use vec_utils
  use UTILS_IO, only: GETUNIT
  use prec_hire
  use saxs_defs
  implicit none
  ! A (somewhat) protected UNIT for file I/O
  integer :: QFILE

  ! SAXS Curve Parameters
  real(kind = real64) :: max_q = 1.0d0
  real(kind = real64) :: saxs_max = 1.0d0
  integer :: num_points = 200
  real(kind = real64) :: delta_q
  real(kind = real64) :: delta_w

  ! Explicit Hydration Parameters
  real(kind = real64) :: dx
  integer :: n_shells

  ! Use and Definition of Target Curve
  logical :: use_target = .true.
  real(kind = real64), dimension(:), allocatable :: target_curve
  real(kind = real64) :: target_curve_relative_threshold  = 0.05 ! Threshold used in check_SAXS_consistency() to detect 
                                                                 ! inconsistency between target_curves and structures

  ! Use and Definition of mean correction Curve
  logical :: use_mean_correction = .true.
  real(kind = real64), dimension(:), allocatable :: mean_correction_curve
  character*100 :: mean_correction_file

  ! Genetic Algorithms Parameters
  real(kind = real64) :: survival_rate ! Proportion of the population that will be discarded
  integer :: population_count = 0
  integer :: SAXS_norm_type = 2 ! 1 for 1-norm, 2 for 2-norm, 3 for infinity-norm

  ! Computation Schemes
  logical :: in_solution_curve = .true. ! .true for In Solution curve, .false. for In Vacuo
  logical :: explicit_sol_contribution = .false. !.true.: explicit removal of solvent contr., .false.: implicit via corrected FF
  logical :: refine_hydration_layer = .false. ! .true. for adding an extra hydration layer with enhanced electron density, .false. otherwise
  logical :: linear_intensity = .false.

  ! FA / CG parameters
  logical :: coarse_grained = .true. ! .true. for CG, .false. for FA (Full-Atom) otherwise
  character*30 :: parameter_file = 'SAXS_grains.dat' ! Either SAXS_grains.dat or SAXS_atoms.dat

  ! Solvent Parameters
  character*4 :: solvent_name = 'HOH '
  real(kind = real64) :: solvent_contrast = 0.0
  real(kind = real64) :: solvent_electron_density = 0.334
  real(kind = real64) :: SAXS_w_shell = 0.003

  ! Module Private Variables
  real(kind = real64) :: max_distance = 100.0d0
  real(kind = real64), dimension(:), allocatable :: sinc_A ! Tabulated cardinal sine function
  integer, parameter :: w_max = 50000
  integer :: MAXSOL = 5000 ! Slightly dynamic, will be increased if necessary
  integer :: num_atoms ! A raw copy of NATOMS (defined in module defs)

  ! Grain specific data
  character*5, dimension(:), allocatable :: Grains
  integer :: grain_number = 0
  real(kind = real64), dimension(:,:), allocatable :: F_grains
  real(kind = real64), dimension(:), allocatable :: Grains_cutoff,Grains_excludedVol,Grains_excludedRadii,Grains_electron_density

  ! Other
  character*1, parameter :: comment_pattern='#'
  logical :: regen_corrected_FF = .false.
  character*180 :: saxs_input_file
  real(kind = real64) :: saxs_alpha = 200.0 ! Coefficient used in Energy Term, in KJ.(...)
  integer :: max_q_point

  contains

    !> @brief Setter method for the parameters of the entire module.
    !> Used for setting up simulation parameters, as well as whecking the existence/generating form factors.
    subroutine set_SAXS_scoring()
      use var_defs, only: nparticles
      implicit none
      integer :: point_q
      logical :: found_mean_correction = .false.
      real(kind = real64) :: q

      !new dummy variable 
      INTEGER :: N_REPLICA = 0

      QFILE=GETUNIT() !find free io unit number
      num_atoms = nparticles
      delta_q = max_q/num_points
      population_count = N_REPLICA
      max_q_point = int(saxs_max / delta_q)

      if (use_target) then
        ! Allocate target SAXS with target file ...
        allocate(target_curve(0:num_points-1))
        ! ... and fill it
        open(unit=QFILE, file='saxs_target.dat', status='old')
        q_loop: do point_q=0, num_points -1
          read(QFILE,*) q, target_curve(point_q)
        end do q_loop
        close(QFILE)
      end if

      select_correction_file: if (in_solution_curve) then
        mean_correction_file = "saxs_mean_correction_Solution.dat"
      else
        mean_correction_file = "saxs_mean_correction_Vacuo.dat"
      end if select_correction_file
      ! Then, ensure its existence
      inquire(file=mean_correction_file, exist=found_mean_correction)
      ! Abort if asked for the correction but can't find the correcting file
      if (use_mean_correction .and. .not. found_mean_correction) then
        print *, 'Could not find mean correction file ', trim(mean_correction_file)
        STOP 5
      end if

      use_mean_correction = use_mean_correction .and. found_mean_correction .and. coarse_grained ! Correction only for CG curves
      if (use_mean_correction ) then
        ! Allocate mean correction SAXS with mean correction file ...
        allocate(mean_correction_curve(0:num_points-1))
        ! ... and fill it
        open(unit=QFILE, file=mean_correction_file, status='old')
        do point_q=0, num_points -1
          read(QFILE,*) q, mean_correction_curve(point_q)
        end do
        close(QFILE)
      end if

      ! Depending on calculation type (FA, CG), change the input parameter file
      if (coarse_grained) then
        parameter_file = 'SAXS_grains.dat'
        call init_hash_CG()
      else
        parameter_file = 'SAXS_atoms.dat'
        call init_hash_FA()
      end if

      ! Init hashes and arrays
      !call tabulate_sinc()

    end subroutine


    !> @brief This method, which should only be called once, ensures that the structures on which
    !> we will be calculating SAXS profiles have a number of electrons matching the target ones
    !> (otherwise, there is obviously an error either in the target or in the structure, but trying
    !> to match them makes no sense and should be avoided)
    subroutine check_SAXS_consistency(pos, ATOMIC_TYPE)
      real(kind = real64), dimension(1:3*num_atoms), intent(in) :: pos
      character*5, dimension(1:num_atoms), intent(in) :: ATOMIC_TYPE
      real(kind = real64), dimension(0:max_q_point-1) :: curve1
      real(kind = real64) :: Esaxs, F_saxs(3*num_atoms)

      call fct_generate_SAXS_curve(pos, ATOMIC_TYPE, curve1, Esaxs, F_saxs)
      !curve1(0:num_points-1) = fct_generate_SAXS_curve(pos, ATOMIC_TYPE) ! ATOMIC_TYPE is sourced in 'module defs'
      if ( abs(10**curve1(0) - 10**target_curve(0) ) >= target_curve_relative_threshold * 10**target_curve(0) ) then
        print '(A,ES15.7,A,ES15.7,A)', 'check_SAXS_consistency(): The target SAXS profile (',10**target_curve(0), &
         & ') does not match the calculated SAXS profile (', 10**curve1(0), &
         & ') in q=0 which is symptomatic of something being very wrong. SAXS_scoring will abort the calculations.'
        STOP 5
      end if
    end subroutine check_SAXS_consistency

    !> @brief Called by saxs_header to output SAXS options to FSAXS (saxs_convergence.dat)
    subroutine write_SAXS_options_to_unit(FUNIT)
    implicit none
    integer, intent(in) :: FUNIT

    write(FUNIT,'(A35,I12   )') 'Number of SAXS Grains found  :     ', grain_number
    write(FUNIT,'(A35,A12   )') 'Curve Type                   :     ', merge(' In Solution', '    In Vacuo', in_solution_curve)
    if ( in_solution_curve .and. refine_hydration_layer ) then
      write(FUNIT,'(A35,f13.10)') 'Dummy Water Radius           :     ', dx
      write(FUNIT,'(A35,I12   )') 'Dummy Water Shells           :     ', n_shells
    end if

    end subroutine write_SAXS_options_to_unit

    !> @brief Call this method to update the SAXS score of a conformation.
!    subroutine compute_score(conformation, write_to_unit)
!      implicit none
      !TODO: commented out to allow compilation
      ! Stand-alone HiRE has no conformations saved
      !type(t_conformations), intent(inout) :: conformation
!      integer, optional :: write_to_unit
!      real(kind = real64), dimension(0:max_q_point-1) :: curve1
!      real(kind = real64) :: Esaxs
!      real(kind = real64), dimension(1:3*num_atoms) :: F_saxs

      ! Generate SAXS Curve
      !TODO: commented out to allow compilation
      ! Stand-alone HiRE has no conformations saved
      !call fct_generate_SAXS_curve(conformation%pos, ATOMIC_TYPE, curve1, Esaxs, F_saxs)
      !curve1(0:num_points-1) = fct_generate_SAXS_curve(conformation%pos, ATOMIC_TYPE) ! ATOMIC_TYPE is sourced in 'module defs'

      ! Compute the new score, and set it
      !TODO: commented out to allow compilation
      ! Stand-alone HiRE has no conformations saved
      !conformation%score = norm_curves(curve1, target_curve)
!
!      if ( present(write_to_unit)) call write_SAXS_curve_to_unit(curve1,write_to_unit)
!
!    end subroutine compute_score


    !> @brief Returns the (log10) SAXS intensity of a given structure (pos + atomNames).
    !function fct_generate_SAXS_curve(pos, ATOMIC_TYPE) result(I_tot)
    subroutine fct_generate_SAXS_curve(pos, ATOMIC_TYPE, logI, Esaxs, F_saxs)

      implicit none
      real(kind = real64), dimension(1:3*num_atoms) :: pos
      character*5, dimension(1:num_atoms) :: ATOMIC_TYPE
      character*5, dimension(:), allocatable :: AtomNames_TOT
      real(kind = real64), dimension(0:max_q_point-1) :: logI, I1, I0
      real(kind = real64) :: max_val
      real(kind = real64), dimension(1:3*num_atoms) :: F_saxs
      real(kind = real64) :: r(3), F_grain(3), r2, qr, qphys, cscale, Esaxs, E_den, E_num, cscale_num, cscale_den, F_pre
      integer :: i, j
      real(kind = real64), dimension(:), allocatable :: pos_SOL, pos_TOT ! TOT = SYSTEM + SOL
      real(kind = real64), dimension(:,:), allocatable :: DistanceMatrix_TOT, F_CG_TOT
      integer :: grain_i, grain_j, q, num_TOT, num_SOL
      real(kind = real64) :: time1, time2
      
      call cpu_time(time1)


      if_refine_hydration_layer: if (in_solution_curve .and. refine_hydration_layer) then
        allocate( pos_SOL(1:3*MAXSOL) )
        pos_SOL = 0.0 ! SAXS-On-The-Fly : Is the initialization necessary ?

        ! Hydrate in pos_SOL
        call hydrate_pos_outer_shell(pos, ATOMIC_TYPE, num_SOL, pos_SOL)
        num_TOT = num_atoms + num_SOL

        ! Allocate pos_TOT with the new solvent molecules
        allocate( pos_TOT(1:3*num_TOT) )
        pos_TOT(1:3*num_atoms) = pos(1:3*num_atoms)
        pos_TOT(3*num_atoms+1:3*num_TOT) = pos_SOL(1:3*num_SOL)

        ! Fill AtomNames_TOT
        allocate(AtomNames_TOT(num_TOT))
        AtomNames_TOT(1:num_atoms) = ATOMIC_TYPE(1:num_atoms)
        AtomNames_TOT(num_atoms+1:num_TOT) = 'HOH  '

        ! Allocate F_CG_TOT and fill it
        allocate( F_CG_TOT(0:num_points-1,1:num_TOT) )
        call fill_structureFactor_array(num_TOT, AtomNames_TOT, F_CG_TOT )! SAXS-On-The-Fly : Need to resize (and fill) ATOMIC_TYPE

        ! Allocate DistanceMatrix_TOT and fill it
        allocate( DistanceMatrix_TOT(1:num_TOT, 1:num_TOT) )
        call fill_half_distance_matrix(num_TOT,pos_TOT,DistanceMatrix_TOT)

      else ! (.not refine_hydration_layer .or. .not. in_solution_curve)
        num_TOT=num_atoms
        ! Allocate F_CG_TOT and fill it
        allocate( F_CG_TOT(0:num_points-1,1:num_TOT) )
        call fill_structureFactor_array(num_TOT, ATOMIC_TYPE, F_CG_TOT )! SAXS-On-The-Fly : Need to resize (and fill) ATOMIC_TYPE
        ! Allocate DistanceMatrix_TOT and fill it
        allocate( DistanceMatrix_TOT(1:num_TOT, 1:num_TOT) )
        allocate( pos_TOT(1:3*num_atoms) )
        call fill_half_distance_matrix(num_TOT, pos, DistanceMatrix_TOT)
        pos_TOT(1:3*num_atoms) = pos(1:3*num_atoms)

      end if if_refine_hydration_layer

! Initialize vars
      I1 = 0.0d0
      F_saxs = 0.0d0
      cscale = 0.0d0
      cscale_num = 0.0d0
      cscale_den = 0.0d0
      E_den = 0.0d0
      E_num = 0.0d0

! q_loop_0
      do grain_i = 1, num_TOT
        I1(0) = I1(0) + F_CG_TOT(0,grain_i)**2
        do grain_j = grain_i+1, num_TOT
          I1(0) = I1(0) + 2 * F_CG_TOT(0,grain_i) * F_CG_TOT(0,grain_j)
        end do
      end do

      q_loop: do q = 1, max_q_point - 1
        do grain_i = 1, num_TOT
! Add diagonal term (i,i) once
          I1(q)  = I1(q) + F_CG_TOT(q,grain_i)**2
! Add half non-diagonal terms twice
          do grain_j = grain_i+1, num_TOT
             I1(q)  = I1(q) + 2 * F_CG_TOT(q,grain_i) * F_CG_TOT(q,grain_j) * &
                 & sin(DistanceMatrix_TOT(grain_j,grain_i) * q * delta_q) / &
                 & (q * delta_q * DistanceMatrix_TOT(grain_j,grain_i)) 
          end do
        end do
      end do q_loop

      ! Log10 everything
      logI(:) = log10(I1(:))
      ! Define linear TARGET curve
      I0(:) = 10**target_curve(:) 
      ! If there is a correction to use, do it
      if (use_mean_correction) then
        logI(:) = logI(:) + mean_correction_curve(:)
        I1(:) = 10**logI(:)
      end if

      ! Compute E_saxs and F_saxs            
      cscale = 1
      do q = 1, max_q_point - 1
        qphys = q * delta_q
        E_num = E_num + ((cscale * I1(q) - I0(q)) * qphys)**2
        F_pre = - 2 * qphys**2 * (cscale * I1(q) - I0(q)) / (max_q_point * I1(0)**2) ! NORMALIZATION ?
      
         do grain_i = 1, num_atoms
            do grain_j = grain_i+1, num_tot
               r = pos_TOT(grain_i*3-2:grain_i*3) - pos_TOT(grain_j*3-2:grain_j*3)
               r2 = dot_product(r,r)
               qr = qphys * dsqrt(r2)
               F_grain = ( r / r2 ) * F_CG_TOT(q, grain_i) * F_CG_TOT(q, grain_j) * ( cos(qr) - sin(qr) / qr )
               F_saxs(grain_i*3-2:grain_i*3) = F_saxs(grain_i*3-2:grain_i*3) + 2 * F_pre * F_grain
               if (grain_j .le. num_atoms) then
                   F_saxs(grain_j*3-2:grain_j*3) = F_saxs(grain_j*3-2:grain_j*3) - 2 * F_pre * F_grain
               endif
            enddo
         enddo
      enddo
      
      Esaxs = E_num / (max_q_point * I1(0)**2)  ! NORMALIZATION ?
        
      ! Deallocate everything
      if (in_solution_curve .and. refine_hydration_layer) then
        deallocate(F_CG_TOT, AtomNames_TOT, pos_SOL, DistanceMatrix_TOT, pos_TOT)
      else
        deallocate(F_CG_TOT, DistanceMatrix_TOT)
      end if


      call cpu_time(time2)
      print *, 'fct_generate_SAXS_curve : ',time2-time1

      return

    !end function fct_generate_SAXS_curve
    end subroutine fct_generate_SAXS_curve

    subroutine write_SAXS_curve_to_file(curve1, file_name)
      implicit none
      real(kind = real64), dimension(0:max_q_point-1), intent(in) :: curve1
      character*360, intent(in) :: file_name
      integer :: q_point,  curvef

      !print *, 'file_name : ', file_name
      curvef = getunit()
      open(unit=curvef,file=file_name,action='write',position='append')

      do q_point = 0, max_q_point-1
       if (linear_intensity) then
         write (curvef,'(f8.3,ES15.7)') q_point*delta_q, 10**curve1(q_point)
       else
         write (curvef,'(f8.3,ES15.7)') q_point*delta_q, curve1(q_point)
       end if
      end do
      close(curvef)

    end subroutine write_SAXS_curve_to_file

    subroutine write_SAXS_curve_to_unit(curve1, unit_number)
      implicit none
      real(kind = real64), dimension(0:max_q_point-1), intent(in) :: curve1
      integer, intent(in) :: unit_number
      integer :: q_point

      do q_point = 0, max_q_point-1
       if (linear_intensity) then
         write (unit_number,'(f8.3,ES15.7)') q_point*delta_q, 10**curve1(q_point)
       else
         write (unit_number,'(f8.3,ES15.7)') q_point*delta_q, curve1(q_point)
       end if
      end do
      write (unit_number,*) ' '

    end subroutine write_SAXS_curve_to_unit

    ! Compute the appropriate norm using SAXS_norm_type
    function norm_curves (curve1, curve2) result(norm)
      implicit none
      real(kind = real64), dimension(0:max_q_point-1) :: curve1, curve2
      real(kind = real64) :: norm
      select case (SAXS_norm_type)
        case (1)
          norm = norm1_curves(curve1, curve2)
        case (2)
          norm = norm2_curves(curve1, curve2)
        case (3)
          norm = norm_R_curves(curve1, curve2)
        case default
          STOP 'Undefined norm type. Aborting.'
          norm = 0.0d0
      end select
    end function norm_curves


    !> @brief Straightforward, norme-2 distance computation between two curves.
    function norm2_curves( curve1, curve2 ) result(norm)
      implicit none
      integer :: q
      real(kind = real64), dimension(0:max_q_point-1) :: curve1, curve2
      real(kind = real64) :: norm

      norm = 0.0d0
      do q = 1, max_q_point - 1
        norm = norm + ((10**curve1(q) - 10**curve2(q)) / (q * delta_q))**2
      end do
      norm = norm / ( max_q_point * (10**curve1(0))**2 ) ! NORMALIZATION
      return
    end function norm2_curves

    !> @brief Straightforward, norme-1 distance computation between two curves.
    function norm1_curves( curve1, curve2 ) result(norm)
      implicit none
      integer :: q
      real(kind = real64), dimension(0:max_q_point-1) :: curve1, curve2
      real(kind = real64) :: norm

      norm = 0.0d0
      do q = 0, max_q_point - 1
        norm = norm + dabs( curve1(q) - curve2(q) )
      end do

      norm = norm / max_q_point ! NORMALIZATION
      return
    end function norm1_curves

    !> @brief Straightforward, infinity-norme distance computation between two curves.
    function norm_R_curves( curve1, curve2 ) result(norm)
      implicit none
      integer :: q
      real(kind = real64), dimension(0:max_q_point-1) :: curve1, curve2
      real(kind = real64) :: norm, norm_num, norm_den, sc, sc_num, sc_den

      norm_num = 0.0d0
      norm_den = 0.0d0
      sc_num = 0.0d0
      sc_den = 0.0d0

     !do q = 1, max_q_point - 1
     !  sc_num = sc_num + ( q * delta_q )**2 * 10**curve1(q) * 10**curve2(q) 
     !  sc_den = sc_den + ( q * delta_q )**2 * ( 10**curve1(q) )**2 
     !end do
     !sc = sc_num / sc_den
      sc = 1 ! Watch out for scale factor

      do q = 1, max_q_point - 1
        norm_num = norm_num + ( q * delta_q * ( sc * 10**curve1(q) - 10**curve2(q) ) )**2
        norm_den = norm_den + ( q * delta_q )**2 * sc * 10**curve1(q) * 10**curve2(q) 
      end do

      norm = norm_num / norm_den 

      return
    end function norm_R_curves


    !> @brief Used to combine two set of coordinates to produce a third one
    !> by linear, symetric mixing of every coordinates.
    subroutine mix_coords(XCoords1, XCoords2, XCoords3, increment_pop_count)
      implicit none
      logical, intent(in) :: increment_pop_count
      real(kind = real64), dimension(1:3*num_atoms), intent(in) :: XCoords1, XCoords2
      real(kind = real64), dimension(1:3*num_atoms), intent(out) :: XCoords3

      XCoords3 = 0.5d0 * (XCoords1 + XCoords2)

      if (increment_pop_count) population_count = population_count + 1
    end subroutine mix_coords


    !> @brief Same as mix_coords, but is a function instead of being a subroutine.
    function fct_mix_coords(XCoords1, XCoords2, increment_pop_count) result(XCoords3)
      implicit none
      logical :: increment_pop_count
      real(kind = real64), dimension(1:3*num_atoms) :: XCoords1, XCoords2
      real(kind = real64), dimension(1:3*num_atoms) :: XCoords3

      XCoords3 = 0.5d0 * (XCoords1 + XCoords2)

      if (increment_pop_count) population_count = population_count + 1

      return
    end function fct_mix_coords


    !> @brief Quicksort algorithm on conf%score.
    !> Taken from http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
!    recursive subroutine QSort_score(a,na)


      ! DUMMY ARGUMENTS
!      integer, intent(in) :: nA
!      type (t_conformations), dimension(nA), intent(in out) :: A

      ! LOCAL VARIABLES
!      integer :: left, right
!      real(kind = real64) :: random
!      real(kind = real64) :: pivot
!      type (t_conformations) :: temp
!      integer :: marker

!          if (nA > 1) then

      !TODO: commented out to allow compilation
      ! Stand-alone HiRE has no conformations saved
!              call random_number(random)
!              pivot = A(int(random*dble(nA-1))+1)%score   ! random pivor (not best performance, but avoids worst-case)
!              left = 0
!              right = nA + 1
!
!              do while (left < right)
!                  right = right - 1
!                  do while (A(right)%score > pivot)
!                      right = right - 1
!                  end do
!                  left = left + 1
!                  do while (A(left)%score < pivot)
!                      left = left + 1
!                  end do
!                  if (left < right) then
!                    temp = A(left)
!                    A(left) = A(right)
!                    A(right) = temp
!                  end if
!              end do
!
!              if (left == right) then
!                  marker = left + 1
!              else
!                  marker = left
!              end if
!
!              call QSort_score(A(1:marker-1),marker-1)
!              call QSort_score(A(marker:nA),nA-marker+1)
!
!          end if
!
!    end subroutine QSort_score


    !> @brief Bubble algorithm on conf%score.
    !> Taken from http://rosettacode.org/wiki/Sorting_algorithms/Bubble_sort#Fortran
!    subroutine BubbleSort_score(A, nA)


      ! DUMMY ARGUMENTS
!      integer, intent(in) :: nA
!      type(t_conformations), dimension(nA), intent(in out) :: A

!      INTEGER :: i, j
!      LOGICAL :: swapped = .TRUE.
!      type(t_conformations) :: temp

      !TODO: commented out to allow compilation
      ! Stand-alone HiRE has no conformations saved
!      DO j = nA-1, 1, -1
!        swapped = .FALSE.
!        DO i = 1, j
!          IF (a(i)%score > a(i+1)%score) THEN
!            temp = a(i)
!            a(i) = a(i+1)
!            a(i+1) = temp
!            swapped = .TRUE.
!          END IF
!        END DO
!        IF (.NOT. swapped) EXIT
!      END DO
!    end subroutine BubbleSort_score

    !> @brief (Re)generates corrected Form Factors for grains.
    !> @warning When multiple threads are launch on the same machine, there seem to be
    !> overwriting problems between them rendering the files corrupted.
    !> To be looked into.
    subroutine gen_corrected_CG_structureFactors()
      use hash_int
      implicit none

      integer :: grain_i, point, solvent_grain_hash
      integer :: unit1, unit2
      logical :: file_exists
      character*100 :: corrected_dat_file, uncorrected_dat_file
      real(kind = real64) :: V_i, q, F_q

      call hash_get(solvent_name, solvent_grain_hash)

      unit1=getunit()
      unit2=getunit()

      regen_corrected_CG_FF: do grain_i=1, grain_number
        ! For each grain ...
        corrected_dat_file ='dat/ff' // trim(Grains(grain_i)) // '.cor.awk.dat'
        inquire(file=corrected_dat_file, exist=file_exists)
        regen: if ( (.not. file_exists) .or. regen_corrected_FF ) then
          ! If the corrected file doesn't exist, or if we force regeneration ...
          uncorrected_dat_file = 'dat/ff' // trim(Grains(grain_i)) // '.awk.dat'
          open(unit=unit1, file=corrected_dat_file, form='formatted')
          is_solvent: if (Grains_electron_density(grain_i) /= 0.0d0) then
            is_not_outer_shell: if (Grains(grain_i)(1:4) .eq. solvent_name) then
              open(unit=unit2, file=uncorrected_dat_file, status='old')
              ! Don't correct the solvent !
              q_loop_1: do point=0, num_points -1
                read(unit2,'(F8.3,F11.8)') q, F_q
                write(unit1,'(F8.3,F13.8)') q, F_q
              end do q_loop_1
            else ! Is outer shell
              uncorrected_dat_file =  'dat/ff' // trim(solvent_name) // '.awk.dat'
              open(unit=unit2, file=uncorrected_dat_file, status='old')
              ! Scale up the outer shell FF with solvent_contrast
              q_loop_2: do point=0, num_points -1
                read(unit2,'(F8.3,F11.8)') q, F_q
                write(unit1,'(F8.3,F13.8)') q, F_q * solvent_contrast
              end do q_loop_2
              uncorrected_dat_file = 'dat/ff' // trim(Grains(grain_i)) // '.awk.dat'
              print *, 'Running : ', 'cp ' // trim(corrected_dat_file) // ' ' // trim(uncorrected_dat_file)
              call system('cp ' // trim(corrected_dat_file) // ' ' // trim(uncorrected_dat_file) )
            end if is_not_outer_shell
          else ! Is not solvent
            open(unit=unit2, file=uncorrected_dat_file, status='old')
            !V_i = (4.0/3.0) * 3.14159 * Grains_Radii(grain_i)**3
            V_i = Grains_ExcludedVol(grain_i)
            q_loop_3: do point=0, num_points -1
              read(unit2,'(F8.3,F11.8)') q, F_q
              write(unit1,'(F8.3,F13.8)') q, F_q - V_i * Grains_electron_density(solvent_grain_hash) * &
                exp( -V_i**(2.0/3) * q**2 / (4.0 * 3.141592653589) )
            end do q_loop_3
          end if is_solvent
          close(unit1)
          close(unit2)
        end if regen
      end do regen_corrected_CG_FF

    end subroutine gen_corrected_CG_structureFactors


    !> @brief Used to speed-up calculation by storing distances between all atoms (time-memory trade-off).
    subroutine fill_half_distance_matrix(num_grain,pos,DistanceMatrix)

      implicit none
      integer, intent(in) :: num_grain
      integer :: i,j
      real(kind = real64), dimension(1:num_grain,1:num_grain), intent(out) :: DistanceMatrix
      real(kind = real64), dimension(1:3*num_grain) ::  pos
      real(kind = real64) :: r_ij

      DistanceMatrix=0.d0
      i_loop: do i=1,num_grain
        j_loop: do j=i+1,num_grain
          r_ij = sqrt((pos(3*i-2)-pos(3*j-2))**2+(pos(3*i-1)-pos(3*j-1))**2+(pos(3*i)-pos(3*j))**2)
          DistanceMatrix(i,j) =r_ij
          DistanceMatrix(j,i)=r_ij     
        end do j_loop
      end do i_loop

    end subroutine fill_half_distance_matrix

    subroutine getLineNumFromUnit_RemoveComments(nf, line_number)
      implicit none

      integer, intent(in) :: nf
      integer,intent(out) :: line_number

      character*100 line
      integer :: stat

      line_number = 0
      do
        read(nf,'(a)',iostat=stat) line
        if (stat /= 0) exit
        if ( line(1:1) /= comment_pattern ) then
          line_number = line_number + 1
        end if
      end do

      return
    end subroutine getLineNumFromUnit_RemoveComments

    !> @brief Main hydration routine. Interesting parameters (dx, n_shells) are module variables that can be
    !> modified with care (being careful of the dx**3 dependance of the algorithm ...)
    !> @details We rely heavily on arrays to speed up computations.
    subroutine hydrate_pos_outer_shell(XCoords, AtomNames, num_SOL, Xcoords_SOL)

      implicit none

      real(kind = real64), intent(in), dimension(1:3*num_atoms) :: XCoords
      character*5, intent(in), dimension(1:num_atoms) :: AtomNames
      real(kind = real64), dimension(1:num_atoms) :: CutoffArrayInternal_pow2, CutoffArrayExternal_pow2
      integer, intent(out) :: num_SOL
      real(kind = real64), intent(out), dimension(1:3*MAXSOL) :: XCoords_SOL

      integer :: lattice_x,lattice_y,lattice_z

      integer :: atom_i

      real(kind = real64) :: r_ij
      real(kind = real64) :: x,y,z
      real(kind = real64) :: cut_off_internal_2, cut_off_external_2
      real(kind = real64), dimension(1:6) :: MinMaxCoordsArray
      integer, dimension(1:6) :: latticeArray

      integer :: contact_num, outunit

      ! Fill the cutoff array with the value of each grain
      call fill_cutoff_array(AtomNames, CutoffArrayInternal_pow2, CutoffArrayExternal_pow2)

      ! Compute dimension of the box ! Format of the array : (/minX,minZ,minZ,maxX,maxY,maxZ/)
      call box_dimensions(XCoords, MinMaxCoordsArray, sqrt( maxval(CutoffArrayExternal_pow2,num_atoms))  ) 

      !Compute max dimensions for the x,y,z loops
      latticeArray = int(MinMaxCoordsArray/dx)

      ! (Re)Initialize num_SOL
      num_SOL = 0

      x_loop: do lattice_x=latticeArray(1), latticeArray(2)
        x = lattice_x*dx
        y_loop: do lattice_y=latticeArray(3), latticeArray(4)
          y = lattice_y*dx
          z_loop: do lattice_z=latticeArray(5), latticeArray(6)
            z = lattice_z*dx
            contact_num = 0
            ! Ensure that there are no atoms within a given cut-off radius
            atom_loop: do atom_i=1,num_atoms
              cut_off_internal_2 = CutoffArrayInternal_pow2(atom_i)
              cut_off_external_2 = CutoffArrayExternal_pow2(atom_i)
              r_ij = (XCoords(3*atom_i-2)-x)**2
              if ( r_ij <= cut_off_external_2 ) then
                r_ij = r_ij + (XCoords(3*atom_i-1)-y)**2
                if ( r_ij <= cut_off_external_2 ) then
                  r_ij= r_ij + (XCoords(3*atom_i)-z)**2
                  ! If clash, abort adding procedure
                  clash: if(r_ij <= cut_off_internal_2 ) then
                    contact_num = 0
                    exit
                  end if clash
                  contact: if (r_ij <= cut_off_external_2 ) then
                    contact_num = contact_num + 1
                  end if contact
                end if
              end if
            end do atom_loop
            add: if (contact_num /= 0 ) then
               num_SOL = num_SOL + 1
               XCoords_SOL(3*num_SOL -2) = x
               XCoords_SOL(3*num_SOL -1) = y
               XCoords_SOL(3*num_SOL   ) = z
            end if add
          end do z_loop
        end do y_loop
      end do x_loop

    end subroutine hydrate_pos_outer_shell

    ! Format of the array : (/minX,minZ,minZ,maxX,maxY,maxZ/)
    subroutine box_dimensions(XCoords, MinMaxCoordsArray, padding_length)

      implicit none
      real(kind = real64),intent(in) :: padding_length ! In nm
      integer :: atom_i
      real(kind = real64),dimension(1:6),intent(out) ::  MinMaxCoordsArray
      real(kind = real64)  XCoords(3*num_atoms)
      integer :: i

      ! lm759 Straight and simple fortran functions to make sure it works
      MinMaxCoordsArray=(/ Minval(XCoords(::3)), Maxval(XCoords(::3)), &
      & Minval(XCoords(2::3)), Maxval(XCoords(2::3)), &
      & Minval(XCoords(3::3)), Maxval(XCoords(3::3)) /)

      ! Init the array 
      ! PDB is limited to 3 digits before the dot, so [-1000..1000] is the range
      MinMaxCoordsArray=(/1000.0,-1000.0,1000.0,-1000.0,1000.0,-1000.0/)

      ! Init the array 
      ! PDB is limited to 3 digits before the dot, so [-1000..1000] is the range
     ! MinMaxCoordsArray=(/1000.0,-1000.0,1000.0,-1000.0,1000.0,-1000.0/)
     ! do atom_i=1,num_atoms ! This algoritm DOES NOT always WORK ! lm759
     !   ! Unrolling manually for performance
     !   if ( XCoords(atom_i*3-2) <= MinMaxCoordsArray(1) ) then   !X_min
     !     MinMaxCoordsArray(1) = XCoords(atom_i*3-2)
     !   else if ( XCoords(atom_i*3-2) >= MinMaxCoordsArray(2) ) then !X_max
     !     MinMaxCoordsArray(2) = XCoords(atom_i*3-2)
     !   end if
     !   if ( XCoords(atom_i*3-1) <= MinMaxCoordsArray(3) ) then   !Y_min
     !     MinMaxCoordsArray(3) = XCoords(atom_i*3-1)
     !   else if ( XCoords(atom_i*3-1) >= MinMaxCoordsArray(4) ) then !Y_max
     !     MinMaxCoordsArray(4) = XCoords(atom_i*3-1)
     !   end if
     !   if ( XCoords(atom_i*3) <= MinMaxCoordsArray(5) ) then   !Z_min
     !     MinMaxCoordsArray(5) = XCoords(atom_i*3)
     !   else if ( XCoords(atom_i*3) >= MinMaxCoordsArray(6) ) then !Z_max
     !     MinMaxCoordsArray(6) = XCoords(atom_i*3)
     !   end if
     ! end do

      ! Pad by padding_length nm in every direction, just to be safe
      ! Warning : this value should be modified for really big coarse grained solvent molecule
      do i=1,6
        ! Substract padding_length to mins, adds it to maxs
        MinMaxCoordsArray(i) = merge(MinMaxCoordsArray(i)-padding_length, MinMaxCoordsArray(i)+padding_length, &
        & modulo(i,2) == 1)
      end do

      return

    end subroutine box_dimensions


    !> @brief Initializes the grain array with Form factors : F_grains
    subroutine init_hash_CG()
      use hash_int
      implicit none
      integer :: grain_i, q
      real(kind = real64) :: bufferValue, bufferValue2, bufferValue3, bufferValue4, temp
      character*180 :: datFile, parameter_line
      character*5 :: bufferKey
      integer :: NF
      integer :: grain_hash_num, stat

      NF=getunit()

      !Get grain number from parameter_file (either SAXS_grains.dat or SAXS_atoms.dat)
      open(unit=NF, file=parameter_file, status="old")
      get_grain_number: do
        read(NF,'(a)', iostat=stat) parameter_line
        if (stat /= 0) exit
        if (parameter_line(1:1) /= comment_pattern) grain_number = grain_number + 1
      end do get_grain_number
      close(NF)

      allocate(Grains(1:grain_number))
      allocate(F_grains(0:num_points-1, 1:grain_number))
      allocate(Grains_cutoff(1:grain_number), Grains_excludedRadii(1:grain_number), Grains_excludedVol(1:grain_number), &
        Grains_electron_density(1:grain_number) )

      ! Make the hash with the grain names from parameter_file (either SAXS_grains.dat or SAXS_atoms.dat)
      call hash_init()
      open(unit=NF, file=parameter_file, status="old")
      grain_i = 1
      build_grain_hash: do while (grain_i <= grain_number)
        read(nf,'(a)',iostat=stat) parameter_line
        if (parameter_line(1:1) /= comment_pattern) then
          read(parameter_line,'(a5)') bufferKey
          call hash_set(bufferKey, grain_i)
          Grains(grain_i)=bufferKey
          grain_i = grain_i + 1
        end if
      end do build_grain_hash
      close(NF)

      ! Get all keys back
      !print *, Grains

      ! Fill the arrays : Grains_cutoff, Grains_excludedRadii, Grains_excludedVol and Grains_electron_density
      open(unit=NF, file=parameter_file, status="old")
      grain_i = 1
      build_grain_arrays: do while (grain_i <= grain_number)
        read (nf, '(a)') parameter_line
        if (parameter_line(1:1) /= comment_pattern) then
          read(parameter_line,'(a5,f4.1,f4.1,f6.1,f6.1,1x,a1)') &
            bufferKey, bufferValue, bufferValue2, bufferValue3, bufferValue4
          call hash_get(bufferKey, grain_hash_num)
          Grains_excludedRadii(grain_hash_num) = bufferValue
          Grains_cutoff(grain_hash_num) = bufferValue2
          Grains_excludedVol(grain_hash_num) = bufferValue3
          Grains_electron_density(grain_hash_num) = bufferValue4
          get_solvent_contrast: if (bufferKey .eq. solvent_name(1:3) // 'o' ) then
            call hash_get(solvent_name, grain_hash_num)
            if (grain_hash_num .eq. 0 ) STOP 'Must put solvent grain before outer solvent grain in SAXS_grains.dat. Aborting.'
              solvent_contrast = bufferValue4 / Grains_electron_density(grain_hash_num)
          end if get_solvent_contrast
          grain_i = grain_i + 1
        end if
      end do build_grain_arrays
      close(NF)

      ! Generate corrected Form Factors
      !call gen_corrected_CG_structureFactors() ! This is now being taken care of in ../saxs/Scattering/ with a make

      ! Build F_grains from the dat/ff* files
      read_dat_files: do grain_i=1, grain_number
        ! Depending on whether or not we are doing the calculation in solvent, we either take the corrected
        ! or the untouched Form Factors
        if ( in_solution_curve .and. .not. explicit_sol_contribution ) then
          datFile= 'dat/ff' // trim(Grains(grain_i)) // '.cor.awk.dat' ! Solvent-corrected FF
        else
          datFile= 'dat/ff' // trim(Grains(grain_i)) // '.awk.dat' ! In vacuo FF
        end if
        open(unit=NF, file=datFile, status='old')
        copy_q: do q=0,num_points-1
          !print *, 'q:', q
          read(NF,'(F8.3,F10.8)') temp, F_grains(q, grain_i)
        end do copy_q
        close(NF)
        ! Added by lm759 to tune excess electron density of the hydration shell
        if ( trim(Grains(grain_i)) .eq. 'HOH' ) then
          F_grains(:, grain_i) = SAXS_w_shell * F_grains(:, grain_i) 
        endif
      end do read_dat_files

    !print *, Grains_excludedRadii
    !print *, Grains_cutoff
    !print *, Grains_excludedVol

    end subroutine init_hash_CG


    !> @brief Initializes the grain array with Form factors : F_grains
    subroutine init_hash_FA()
      use hash_int
      implicit none
      integer :: grain_i, q_point
      real(kind = real64) :: bufferValue, bufferValue2, V_i, q
      character*180 :: parameter_line
      character*5 :: bufferKey
      integer :: NF
      integer :: grain_hash_num, stat
      ! Cromer-Mann parameters
      real(kind = real64) :: a1, a2, a3, a4, c, b1, b2, b3, b4

      NF=getunit()

      !Get grain number from parameter_file (either SAXS_grains.dat or SAXS_atoms.dat)
      open(unit=NF, file=parameter_file, status="old")
      get_grain_number: do
        read(NF,'(a)', iostat=stat) parameter_line
        if (stat /= 0) exit
        if (parameter_line(1:1) /= comment_pattern) grain_number = grain_number + 1
      end do get_grain_number
      close(NF)

      allocate(Grains(1:grain_number))
      allocate(F_grains(0:num_points-1, 1:grain_number))
      allocate(Grains_cutoff(1:grain_number), Grains_excludedRadii(1:grain_number), Grains_excludedVol(1:grain_number), &
        Grains_electron_density(1:grain_number) )

      ! Make the hash with the grain names from parameter_file (either SAXS_grains.dat or SAXS_atoms.dat)
      call hash_init()
      open(unit=NF, file=parameter_file, status="old")
      grain_i = 1
      build_grain_hash: do while (grain_i <= grain_number)
        read(nf,'(a)',iostat=stat) parameter_line
        if (parameter_line(1:1) /= comment_pattern) then
          read(parameter_line,'(a5)') bufferKey
          call hash_set(bufferKey, grain_i)
          Grains(grain_i)=bufferKey
          grain_i = grain_i + 1
        end if
      end do build_grain_hash
      close(NF)

      ! DEBUG
      !print *, Grains
      !print *, grain_number

      ! Fill the arrays : Grains_cutoff, Grains_excludedRadii, Grains_excludedVol and Grains_electron_density
      open(unit=NF, file=parameter_file, status="old")
      grain_i = 1
      build_grain_arrays: do while (grain_i <= grain_number)
        read (nf, '(a)') parameter_line
        if (parameter_line(1:1) /= comment_pattern) then
          read(parameter_line,'(a5,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f4.1,f6.2)',iostat=stat) &
            bufferKey, a1, a2, a3, a4, c, b1, b2, b3, b4, bufferValue, bufferValue2
          print '(a5,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f4.1,f6.2)', &
            bufferKey, a1, a2, a3, a4, c, b1, b2, b3, b4, bufferValue, bufferValue2
          call hash_get(bufferKey, grain_hash_num)
          Grains_cutoff(grain_hash_num) = bufferValue
          Grains_excludedVol(grain_hash_num) = bufferValue2
          compute_q: do q_point=0,num_points-1
            q = q_point * delta_q
            implicit_correct_form_factors: if ( in_solution_curve .and. .not. explicit_sol_contribution ) then
              V_i = Grains_excludedVol(grain_hash_num)
              F_grains(q_point, grain_i) = c + a1*EXP(-b1*(q**2)/157.91d0) + a2*EXP(-b2*(q**2)/157.91d0) + &
                                               a3*EXP(-b3*(q**2)/157.91d0) + a4*EXP(-b4*(q**2)/157.91d0) - &
                   V_i * solvent_electron_density * exp( -V_i**(2.0/3) * q**2 / (4.0 * 3.141592653589) )
            else
              F_grains(q_point, grain_i) = c + a1*EXP(-b1*(q**2)/157.91d0) + a2*EXP(-b2*(q**2)/157.91d0) + &
                                               a3*EXP(-b3*(q**2)/157.91d0) + a4*EXP(-b4*(q**2)/157.91d0)
            end if implicit_correct_form_factors
          end do compute_q
          grain_i = grain_i + 1
        end if
      end do build_grain_arrays
      close(NF)

    end subroutine init_hash_FA


    subroutine fill_structureFactor_array(num_grain, GrainName, F_q_CG)
      implicit none
      integer, intent(in) :: num_grain
      character*5, intent(in), dimension(1:num_grain) :: GrainName
      real(kind = real64), intent(out),dimension(0:num_points-1, 1:num_grain ) :: F_q_CG

      if (coarse_grained) then
        call fill_CG_structureFactor_array(num_grain, GrainName, F_q_CG)
      else
        call fill_FA_structureFactor_array(num_grain, GrainName, F_q_CG)
      end if

    end subroutine fill_structureFactor_array

    !> @brief Fills a F_q_CG array with FF values from F_grains
    subroutine fill_CG_structureFactor_array(num_grain, GrainName, F_q_CG)

      use hash_int
      implicit none

      integer, intent(in) :: num_grain
      character*5, intent(in), dimension(1:num_grain) :: GrainName
      real(kind = real64), intent(out),dimension(0:num_points-1, 1:num_grain ) :: F_q_CG

      integer :: grain_hash_num
      integer :: grain_i

      ! Finally, fill the array accordingly
      match_cutoff: do grain_i=1,num_grain
        call hash_get( adjustl(trim(GrainName(grain_i))) , grain_hash_num )
        if (grain_hash_num /= 0 ) then
          F_q_CG(0:num_points-1, grain_i) = F_grains(0:num_points-1, grain_hash_num)
        else
          print *, "Couldn't find structure factor parameters for Grain : "
          print '(i5,a)', grain_i, GrainName(grain_i)
          F_q_CG(0:num_points-1, grain_i) = 0.0
          STOP 5
        end if
      end do match_cutoff

    end subroutine fill_CG_structureFactor_array


    !> @brief Fills a F_q_CG array with FF values from F_grains
    subroutine fill_FA_structureFactor_array(num_grain, GrainName, F_q_CG)

      use hash_int
      implicit none

      integer, intent(in) :: num_grain
      character*5, intent(in), dimension(1:num_grain) :: GrainName
      character*3 :: atom
      real(kind = real64), intent(out),dimension(0:num_points-1, 1:num_grain ) :: F_q_CG

      integer :: grain_hash_num
      integer :: grain_i

      ! Finally, fill the array accordingly
      match_cutoff: do grain_i=1,num_grain
        atom = adjustl(trim(GrainName(grain_i))) ! TODO : Might contain error (for instance, OH2 -> H, not O)
        if (atom(1:1) .eq. 'O') then
          if (atom(3:3) .eq. 'P') then
            ! Account for phosphate's negatively charged oxygen atoms
            atom = 'O-'
          else
            atom = atom(1:1)
          end if
        else
          atom = atom(1:1)
        end if
        call hash_get( atom(1:2) , grain_hash_num )
        if (grain_hash_num /= 0 ) then
          F_q_CG(0:num_points-1, grain_i) = F_grains(0:num_points-1, grain_hash_num)
        else
          print *, "Couldn't find structure factor parameters for Grain : "
          print '(i5,a,a)', grain_i, ' ', atom
          F_q_CG(0:num_points-1, grain_i) = 0.0
          STOP 5
        end if
      end do match_cutoff

    end subroutine fill_FA_structureFactor_array


    subroutine fill_cutoff_array(AtomNames, CutoffArrayInternal_pow2, CutoffArrayExternal_pow2)

      use hash_int
      implicit none

      character*5, intent(in), dimension(1:num_atoms) :: AtomNames
      real(kind = real64), dimension(1:num_atoms), intent(out) :: CutoffArrayInternal_pow2
      real(kind = real64), dimension(1:num_atoms), intent(out) :: CutoffArrayExternal_pow2

      integer :: atom_i, grain_hash_num
      character*5 :: atom_name

      ! Finally, fill the cutoff array accordingly, using * as the unknown value
      match_cutoff: do atom_i=1,num_atoms
          if (coarse_grained) then
            atom_name = trim(AtomNames(atom_i))
          else
            atom_name(1:5) = adjustl(AtomNames(atom_i))
            atom_name(1:5) = atom_name(1:1) // '    '
          end if
          call hash_get( atom_name, grain_hash_num )
          if (grain_hash_num .eq. 0) then
            print *, "Couldn't find parameters for atom/grain : ",atom_name
            STOP 6
          end if
          !if ( grain_hash_num == 0 ) call hash_get('*', grain_hash_num) ! Last chance to get a non-null value
          CutoffArrayInternal_pow2(atom_i) = ( 0.5   * dx                    + Grains_cutoff(grain_hash_num) )**2
          CutoffArrayExternal_pow2(atom_i) = ( 0.99  * (n_shells + 0.5) * dx + Grains_cutoff(grain_hash_num) )**2
      end do match_cutoff

    end subroutine fill_cutoff_array

end module SAXS_scoring
