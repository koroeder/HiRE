MODULE SAXS_CALCS
   USE UTILS_IO, ONLY: GETUNIT
   USE VAR_DEFS, ONLY: NPARTICLES
   USE PREC_HIRE
   USE SAXS_DEFS

   IMPLICIT NONE

   INTEGER :: NQPOINTS = 0                            ! number of q points
   REAL(KIND = REAL64), PARAMETER :: DELTA_Q = 0.002  ! spacing of q points

   !Target curve 
   LOGICAL :: USE_TARGET = .TRUE. ! use target curve
   REAL(KIND = REAL64), ALLOCATABLE :: TARGET_CURVE(:) ! target curve data
   REAL(KIND = REAL64), PARAMETER :: TC_REL_THRESHOLD = 0.05 ! threshold for difference between predicted and target in consistency check
   CHARACTER(LEN=20) :: TC_FILE = "saxs_target.dat"

   !Mean correction curve
   LOGICAL :: USE_MEAN_CORRECTION = .TRUE.  !use mean correction
   REAL(KIND = REAL64), ALLOCATABLE :: MEAN_CORR_CURVE(:) ! correction curve data
   CHARACTER(LEN=40) :: MC_FILE

   !Computational schemes
   LOGICAL :: IN_SOLUTION_CURVE = .TRUE. ! if T, curve is in solution, else in vacuo 

   !Grain data
   CHARACTER(LEN=25) :: PARAMETER_FILE = "SAXS_grains.dat"
   CHARACTER(LEN=1), PARAMETER :: COMMENT_PATTERN = "#"
   INTEGER :: NGRAINS = 0
   CHARACTER(LEN=5), ALLOCATABLE :: GRAINS(:)
   REAL(KIND=REAL64), ALLOCATABLE :: FACTOR_GRAINS(:,:)
   REAL(KIND=REAL64), ALLOCATABLE :: GRAINS_CUTOFF(:)
   REAL(KIND=REAL64), ALLOCATABLE :: GRAINS_EXCLV(:)
   REAL(KIND=REAL64), ALLOCATABLE :: GRAINS_EXCLR(:)
   REAL(KIND=REAL64), ALLOCATABLE :: GRAINS_EDENSITY(:)

   ! Structure factor data for solute
   REAL(KIND=REAL64), ALLOCATABLE :: STRUCT_FACTORS_CG(:,:)
   REAL(KIND=REAL64), ALLOCATABLE :: STRUCT_FACTORS_PRODS_CG(:,:,:)
   CONTAINS

      SUBROUTINE SAXS_SETUP()
         USE UTILS_IO, ONLY: FILE_EXIST
         IMPLICIT NONE
         INTEGER :: QFILE
         INTEGER :: IDX
         REAL(KIND = REAL64) :: DUMMY

         NQPOINTS = INT(SAXSMAX/DELTA_Q)
         !setup SAXS target curve
         IF (USE_TARGET) THEN
            ALLOCATE(TARGET_CURVE(NQPOINTS))
            QFILE = GETUNIT()
            OPEN(UNIT=QFILE, FILE=TC_FILE, STATUS='OLD')
            DO IDX=1,NQPOINTS
               READ(QFILE,*) QDUMMY, TARGET_CURVE(IDX)
            END DO
            CLOSE(QFILE)
         END IF
         !setup mean correction
         ! first set the correct data file and check it is there
         IF (IN_SOLUTION_CURVE) THEN
            MC_FILE = "saxs_mean_correction_Solution.dat"
         ELSE
            MC_FILE = "saxs_mean_correction_Vacuo.dat"
         END IF
         IF (USE_MEAN_CORRECTION.AND.(.NOT.FILE_EXIST(MC_FILE))) THEN
            WRITE(*,*) " saxs_setup> Could not find mean correction file ", MC_FILE
            STOP
         END IF
         !now load data
         IF (USE_MEAN_CORRECTION) THEN
            ALLOCATE(MEAN_CORR_CURVE(NQPOINTS))
            QFILE = GETUNIT()
            OPEN(UNIT=QFILE, FILE=MC_FILE, STATUS='OLD')
            DO IDX =1,NQPOINTS
               READ(QFILE,*) DUMMY, MEAN_CORR_CURVE(IDX)
            END DO
            CLOSE(QFILE)
         END IF
         ! Load the grain data (only CG option)
         CALL INIT_HASH_CG()
         ! TODO: REFACTOR INIT_HASH
         ALLOCATE(STRUCT_FACTORS_CG(NQPOINTS,NPARTICLES),STRUCT_FACTORS_PRODS_CG(NQPOINTS,NPARTICLES,NPARTICLES))
         CALL FILL_STRCTURE_FACTOR_CG()
      END SUBROUTINE SAXS_SETUP


      SUBROUTINE GENERATE_SAXS_CURVE(X, LOGI, ESAXS, FSAXS)
         USE VAR_DEFS, ONLY: ATNAMES => IGRAPH
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NPARTICLES)
         REAL(KIND=REAL64) :: LOGI(NQPOINTS)
         REAL(KIND=REAL64), INTENT(OUT) :: ESAXS
         REAL(KIND=REAL64), INTENT(OUT) :: FSAXS(3*NPARTICLES)

         ! variables for CG particles+solvent
         INTEGER :: NUM_TOT = 0
         REAL(KIND=REAL64), ALLOCATABLE :: X_TOT(:)
         CHARACTER(LEN=4), ALLOCATABLE :: ALLNAMES(:)
         REAL(KIND=REAL64), ALLOCATABLE :: STRUCT_FACTORS(:,:)
         REAL(KIND=REAL64), ALLOCATABLE :: DISTMATRIX(:,:)
         REAL(KIND=REAL64), ALLOCATABLE :: STRUCT_FACTORS_PRODS(:,:,:)
         !solvation variables
         INTEGER :: NUM_SOL = 0
         REAL(KIND=REAL64), ALLOCATABLE :: XSOL(:)

         ! add hydration layer if needed
         IF (IN_SOLUTION_CURVE.AND.REFINE_HYDRATION_LAYER) THEN
            CALL HYDRATE_STRUCTURE(X, ATNAMES, NUM_SOL, XSOL)
            NUM_TOT = NPARTICLES + NUM_SOL
         ELSE
            NUM_TOT = NPARTICLES
         END IF
         ! allocation block for variables
         ALLOCATE(X_TOT(3*NUM_TOT),ALLNAMES(NUM_TOT), &
                  STRUCT_FACTORS(NUM_TOT),DISTMATRIX(NUM_TOT,NUM_TOT), &
                  STRUCT_FACTORS_PRODS(NUM_TOT,NUM_TOT))
         ! copy CG particle data
         X_TOT(1:3*NPARTICLES) = X(1:3*NPARTICLES)
         ALLNAMES(1:NPARTICLES) = ATNAMES(1:NPARTICLES)

         IF (NUM_TOT.GT.NPARTICLES) THEN
            X_TOT(3*NPARTICLES+1:3*NUM_TOT) = XSOL(1:3*NUM_SOL)
            ALLNAMES(NPARTICLES+1:NUM_TOT) = 'HOH'
         END IF
         ! fill structure factor array
         CALL FILL_STRCTURE_FACTOR()
         ! fill distance matrix
         CALL FILL_HALF_DIST_MATRIX()

         ! get products of structure factors
         STRUCT_FACTORS_PRODS(:,:) = 0.0D0
         DO I=1,NUM_TOT
            DO J = I,NUM_TOT
               STRUCT_FACTORS_PRODS(I,J) = STRUCT_FACTORS(I)*STRUCT_FACTORS(J)
               STRUCT_FACTORS_PRODS(J,I) = STRUCT_FACTORS(I,J)
            END DO
         END DO

         I1 = SUM(STRUCT_FACTORS_PRODS)



      END SUBROUTINE GENERATE_SAXS_CURVE

      SUBROUTINE INIT_HASH_CG()
         USE HASH_INT
         IMPLICIT NONE
         INTEGER :: FUNIT, STATUS
         CHARACTER(LEN=150) :: PARAMETER_LINE 

         ! first obtain the number of grains
         FUNIT = GETUNIT()
         OPEN(UNIT=FUNIT, FILE=PARAMETER_FILE, STATUS='OLD')
         DO
            READ(FUNIT, '(A)', IOSTAT=STATUS) PARAMETER_LINE
            IF (STATUS.NE.0) EXIT
            IF (PARAMETER_LINE(1:1).NE.COMMENT_PATTERN) THEN
               NGRAINS = NGRAINS + 1
            END IF
         END DO
         CLOSE(FUNIT)
         ! allocate arrays and initialise hasing
         CALL ALLOCATE_GRAIN_DATA()
         CALL HASH_INIT()

         ! second parse - get grain names hashed
         FUNIT = GETUNIT()
         OPEN(UNIT=FUNIT, FILE=PARAMETER_FILE, STATUS='OLD')



      END SUBROUTINE INIT_HASH_CG

      SUBROUTINE ALLOCATE_GRAIN_DATA()
         IMPLICIT NONE
         CALL DEALLOC_GRAIN_DATA()
         ALLOCATE(FACTOR_GRAINS(NQPOINTS,NGRAINS))
         ALLOCATE(GRAINS(NGRAINS))
         ALLOCATE(GRAINS_CUTOFF(NGRAINS))
         ALLOCATE(GRAINS_EDENSITY(NGRAINS))
         ALLOCATE(GRAINS_EXCLR(NGRAINS))
         ALLOCATE(GRAINS_EXCLV(NGRAINS))
      END SUBROUTINE ALLOCATE_GRAIN_DATA

      SUBROUTINE DEALLOC_GRAIN_DATA()
         IMPLICIT NONE
         IF (ALLOCATED(FACTOR_GRAINS)) DEALLOCATE(FACTOR_GRAINS)
         IF (ALLOCATED(GRAINS)) DEALLOCATE(GRAINS)
         IF (ALLOCATED(GRAINS_CUTOFF)) DEALLOCATE(GRAINS_CUTOFF)
         IF (ALLOCATED(GRAINS_EDENSITY)) DEALLOCATE(GRAINS_EDENSITY)
         IF (ALLOCATED(GRAINS_EXCLR)) DEALLOCATE(GRAINS_EXCLR)
         IF (ALLOCATED(GRAINS_EXCLV)) DEALLOCATE(GRAINS_EXCLV)
      END SUBROUTINE DEALLOC_GRAIN_DATA

END MODULE SAXS_CALCS