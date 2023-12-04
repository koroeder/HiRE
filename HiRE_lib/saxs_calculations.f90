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
      END SUBROUTINE SAXS_SETUP


      SUBROUTINE GENERATE_SAXS_CURVE(X, ATNAMES, LOGI, ESAXS, FSAXS)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NPARTICLES)
         CHARACTER(LEN=4), INTENT(IN) :: ATNAMES(NPARTICLES)
         REAL(KIND=REAL64) :: LOGI(NQPOINTS)
         REAL(KIND=REAL64), INTENT(OUT) :: ESAXS
         REAL(KIND=REAL64), INTENT(OUT) :: FSAXS(3*NPARTICLES)

         ! variables for CG particles+solvent
         INTEGER :: NUM_TOT = 0
         REAL(KIND=REAL64), ALLOCATABLE :: X_TOT(:)
         CHARACTER(LEN=4), ALLOCATABLE :: ALLNAMES(:)
         REAL(KIND=REAL64), ALLOCATABLE :: STRUCT_FACTORS(:)
         REAL(KIND=REAL64), ALLOCATABLE :: DISTMATRIX(:,:)
         REAL(KIND=REAL64), ALLOCATABLE :: STRUCT_FACTORS_PRODS(:,:)
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

END MODULE SAXS_CALCS