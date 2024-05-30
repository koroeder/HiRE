MODULE SAXS_CALCS
   USE UTILS_IO, ONLY: GETUNIT
   USE VAR_DEFS, ONLY: NPARTICLES
   USE PREC_HIRE
   USE SAXS_DEFS

   IMPLICIT NONE

   INTEGER :: NQPOINTS = 0                            ! number of q points
   REAL(KIND = REAL64), PARAMETER :: DELTA_Q = 0.002  ! spacing of q points
   REAL(KIND = REAL64), ALLOCATABLE :: QPOINTS(:)     ! Q point values

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
   LOGICAL :: EXPLICIT_SOL_CONTR = .FALSE. !.true.: explicit removal of solvent contr., .false.: implicit via corrected FF
   LOGICAL :: REFINE_HYDRATION_LAYER = .FALSE. ! .true. for adding an extra hydration layer with enhanced electron density
   LOGICAL :: LINEAR_INTENSITY = .FALSE.

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

   !solvation variables
   INTEGER :: NUM_SOL = 0
   REAL(KIND=REAL64), ALLOCATABLE :: XSOL(:)

   ! Solvent data
   CHARACTER(LEN=5) :: SOLVENT_NAME = 'HOH'
   REAL(KIND=REAL64) :: SOLVENT_CONTRAST = 0.0
   REAL(KIND=REAL64), PARAMETER :: SOLVENT_ELEC_DENSITY = 0.334
   REAL(KIND=REAL64) :: SAXS_W_SHELL = 0.003

   REAL(KIND=REAL64) :: DX
   INTEGER :: N_SHELLS

   CONTAINS

      SUBROUTINE SAXS_SETUP()
         USE UTILS_IO, ONLY: FILE_EXIST
         USE VAR_DEFS, ONLY: ATNAMES => IGRAPH
         IMPLICIT NONE
         INTEGER :: QFILE
         INTEGER :: IDX, I
         REAL(KIND = REAL64) :: DUMMY, QDUMMY

         NQPOINTS = INT(SAXSMAX/DELTA_Q)
         ALLOCATE(QPOINTS(NQPOINTS))
         DO I=1,NQPOINTS
            QPOINTS(I) = (I-1)*DELTA_Q
         END DO

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

         ALLOCATE(STRUCT_FACTORS_CG(NQPOINTS,NPARTICLES),STRUCT_FACTORS_PRODS_CG(NQPOINTS,NPARTICLES,NPARTICLES))
         CALL FILL_STRUCTURE_FACTOR_CG(NPARTICLES,ATNAMES,STRUCT_FACTORS_CG,STRUCT_FACTORS_PRODS_CG)
      END SUBROUTINE SAXS_SETUP


      SUBROUTINE GENERATE_SAXS_CURVE(X, LOGI, ESAXS, FSAXS, GRADT)
         USE VEC_UTILS, ONLY: VEC_DIFF, NORMED_VEC
         USE UTILS_MATHS
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NPARTICLES)
         REAL(KIND=REAL64), INTENT(OUT) :: LOGI(NQPOINTS)
         REAL(KIND=REAL64), INTENT(OUT) :: ESAXS
         REAL(KIND=REAL64), INTENT(OUT) :: FSAXS(3*NPARTICLES)
         LOGICAL, INTENT(IN) :: GRADT

         REAL(KIND=REAL64), ALLOCATABLE :: DISTMATRIX(:,:)
         !variables for calculation of curve
         REAL(KIND=REAL64) :: I1(NQPOINTS), I0(NQPOINTS), I1Q   !intensities
         REAL(KIND=REAL64) :: SF_LOCAL(NPARTICLES,NPARTICLES)   !structure factors for current Q
         REAL(KIND=REAL64) :: ARG, CSCALE, FSCALE, EDUMMY, QPHYS, WC !dummy variables for memorisation
         REAL(KIND=REAL64) :: DIFF, DIFFQ, F_PRE, R2, QR, T1, T2 !dummy variables for memorisation
         REAL(KIND=REAL64) :: R(3), FGRAIN(3), FACTING(3) ! dummy variables for memorisation
         INTEGER :: Q, I, J, IDX, JDX

         !WRITE(*,*) " "
         !WRITE(*,*) "start of routine "
         !call cpu_time(t1)

         ! allocation block for variables
         ALLOCATE(DISTMATRIX(NPARTICLES,NPARTICLES))
         ! fill distance matrix
         CALL FILL_HALF_DIST_MATRIX(NPARTICLES,X,DISTMATRIX)

         !call cpu_time(t2)
         !WRITE(*,*) "After allocations, time taken: ", t2-t1
         !t1 = t2         

         I1(1:NQPOINTS) = 0.0D0
         ! loop for q=0 (note the array index for I1 is 1)
         !DO I = 1,NPARTICLES
         !   I1(1) = I1(1) + STRUCT_FACTORS_PRODS_CG(1,I,I)
         !   DO J = I+1,NPARTICLES
         !      I1(1) = I1(1) + 2*STRUCT_FACTORS_PRODS_CG(1,I,J)
         !   END DO
         !END DO
         I1(1) = SUM(STRUCT_FACTORS_PRODS_CG(1,:,:))

         !call cpu_time(t2)
         !WRITE(*,*) "after qloop_0, time taken: ", t2-t1
         !t1 = t2

         ! loop for q/=0 (starting index is 2 for the I1, QPOINTS and prods arrays)
         DO Q=2,NQPOINTS
            SF_LOCAL(:,:) = STRUCT_FACTORS_PRODS_CG(Q,:,:)
            QPHYS = QPOINTS(Q)
            I1Q = 0.0D0
            DO I= 1, NPARTICLES
               I1Q = I1Q + SF_LOCAL(I,I)
               DO J=I+1,NPARTICLES
                  ARG = DISTMATRIX(J,I) * QPHYS
                  !I1Q = I1Q + 2*SF_LOCAL(I,J)*DSIN(ARG)/ARG
                  I1Q = I1Q + 2*SF_LOCAL(I,J)*MYSIN(ARG)/ARG
               END DO
            END DO
            I1(Q) = I1Q
         END DO

         !call cpu_time(t2)
         !WRITE(*,*) "after qloop, time taken: ", t2-t1
         !t1 = t2

         !apply logarithms
         LOGI(1:NQPOINTS) = DLOG10(I1(1:NQPOINTS))
         !Define target curve
         I0(1:NQPOINTS) = 10**TARGET_CURVE(1:NQPOINTS)
         !Use correction
         IF (USE_MEAN_CORRECTION) THEN
            LOGI(1:NQPOINTS) = LOGI(1:NQPOINTS) + MEAN_CORR_CURVE(1:NQPOINTS)
            I1(1:NQPOINTS) = 10**LOGI(1:NQPOINTS)
         END IF

         !call cpu_time(t2)
         !WRITE(*,*) "after taking log10, time taken: ", t2-t1
         !t1 = t2
         EDUMMY = 0.0D0

         IF (GRADT) THEN
            ! Compute SAXS E and potentially F, we try to use memorisation as much as possible
            CSCALE = I0(1)/I1(1)
            FSCALE = 1.0D0/(NQPOINTS*I1(1)**2)
            DO Q=2,NQPOINTS
               QPHYS = QPOINTS(Q)
               DIFF = CSCALE*I1(Q) - I0(Q)
               !DIFFQ = DIFF*QPHYS
               WC = 0.02 + QPHYS + 5*QPHYS**2
               !EDUMMY = EDUMMY + DIFFQ**2
               EDUMMY = EDUMMY + (DIFF*WC)**2
               !F_PRE = -2*QPHYS*DIFFQ*FSCALE
               F_PRE = -2*FSCALE*DIFF*WC**2
               SF_LOCAL(:,:) = STRUCT_FACTORS_PRODS_CG(Q,:,:)
               DO I = 1,NPARTICLES
                  IDX = 3*I
                  DO J = I+1,NPARTICLES
                     JDX = 3*J
                     R(1:3) = VEC_DIFF(X(IDX-2:IDX),X(JDX-2:JDX))
                     R2 = DOT_PRODUCT(R,R)
                     QR = QPHYS*DSQRT(R2)
                     !FGRAIN = (R/R2)*SF_LOCAL(I,J)*(DCOS(QR)-DSIN(QR)/QR)
                     FGRAIN = (R/R2)*SF_LOCAL(I,J)*(MYCOS(QR)-MYSIN(QR)/QR)
                     FACTING = 2*FGRAIN*F_PRE
                     FSAXS(IDX-2:IDX) = FSAXS(IDX-2:IDX) + FACTING
                     FSAXS(JDX-2:JDX) = FSAXS(JDX-2:JDX) - FACTING
                  END DO
               END DO
            END DO
         ELSE
            CSCALE = I0(1)/I1(1)
            DO Q=2,NQPOINTS
               QPHYS = QPOINTS(Q)
               DIFF = CSCALE*I1(Q) - I0(Q)
               DIFFQ = DIFF*QPHYS
               EDUMMY = EDUMMY + DIFFQ**2
            END DO
            FSAXS(:) = 0.0D0
         END IF
         !call cpu_time(t2)
         !WRITE(*,*) "after force loop, time taken: ", t2-t1
         !t1 = t2

         !normalise
         ESAXS = EDUMMY/(NQPOINTS*I0(1))
         DEALLOCATE(DISTMATRIX)

         WRITE(*,*) "ESAXS", ESAXS
         
         !call cpu_time(t2)
         !WRITE(*,*) "end of routine, time taken: ", t2-t1
         !WRITE(*,*) "--------------------"
         !WRITE(*,*) " "
      END SUBROUTINE GENERATE_SAXS_CURVE

      SUBROUTINE WRITE_SAXS_CURVE_TO_UNIT(CURVE, UNIT)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: CURVE(NQPOINTS)
         INTEGER, INTENT(IN) :: UNIT
         INTEGER :: Q

         DO Q=1,NQPOINTS
            IF (LINEAR_INTENSITY) THEN
               WRITE(UNIT,'(F8.3,ES15.7)') (Q-1)*DELTA_Q, 10**CURVE(Q)
            ELSE
               WRITE(UNIT,'(F8.3,ES15.7)') (Q-1)*DELTA_Q, CURVE(Q)               
            END IF
         END DO
         WRITE(UNIT,'(A)') "  "
      END SUBROUTINE WRITE_SAXS_CURVE_TO_UNIT

      SUBROUTINE GENERATE_SAXS_CURVE_HYDRATE(X, LOGI, ESAXS, FSAXS)
         USE VAR_DEFS, ONLY: ATNAMES => IGRAPH
         USE VEC_UTILS, ONLY: VEC_DIFF, NORMED_VEC
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

         !variables for calculation of curve
         REAL(KIND=REAL64) :: I1(NQPOINTS), I0(NQPOINTS)         !intensities
         REAL(KIND=REAL64) :: ARG, CSCALE, FSCALE, EDUMMY, QPHYS !dummy variables for memorisation
         REAL(KIND=REAL64) :: DIFF, DIFFQ, F_PRE, R2, QR !dummy variables for memorisation
         REAL(KIND=REAL64) :: R(3), FGRAIN(3), FACTING(3) ! dummy variables for memorisation
         INTEGER :: Q, I, J, IDX, JDX

         ! add hydration layer
         CALL HYDRATE_STRUCTURE(X)
         NUM_TOT = NPARTICLES + NUM_SOL

         ! allocation block for variables
         ALLOCATE(X_TOT(3*NUM_TOT),ALLNAMES(NUM_TOT), &
                  STRUCT_FACTORS(NQPOINTS,NUM_TOT),DISTMATRIX(NUM_TOT,NUM_TOT), &
                  STRUCT_FACTORS_PRODS(NQPOINTS,NUM_TOT,NUM_TOT))
         ! copy CG particle data
         X_TOT(1:3*NPARTICLES) = X(1:3*NPARTICLES)
         ALLNAMES(1:NPARTICLES) = ATNAMES(1:NPARTICLES)

         X_TOT(3*NPARTICLES+1:3*NUM_TOT) = XSOL(1:3*NUM_SOL)
         ALLNAMES(NPARTICLES+1:NUM_TOT) = 'HOH'
         CALL FILL_STRUCTURE_FACTOR_CG(NUM_TOT,ALLNAMES,STRUCT_FACTORS,STRUCT_FACTORS_PRODS)

         ! fill distance matrix
         CALL FILL_HALF_DIST_MATRIX(NUM_TOT,X_TOT,DISTMATRIX)

         I1(1:NQPOINTS) = 0.0D0
         ! loop for q=0 (note the array index for I1 is 1)
         !DO I = 1,NUM_TOT
         !   I1(1) = I1(1) + STRUCT_FACTORS_PRODS(1,I,I)
         !   DO J = I+1,NUM_TOT
         !      I1(1) = I1(1) + 2*STRUCT_FACTORS_PRODS(1,I,J)
         !   END DO
         !END DO
         I1(1) = SUM(STRUCT_FACTORS_PRODS(1,:,:))

         ! loop for q/=0 (starting index is 2 for the I1, QPOINTS and prods arrays)
         DO Q=2,NQPOINTS
            DO I= 1, NUM_TOT
               I1(Q) = I1(Q) + STRUCT_FACTORS_PRODS(Q,I,I)
               DO J=I+1,NUM_TOT
                  ARG = DISTMATRIX(J,I) * QPOINTS(Q)
                  I1(Q) = I1(Q) + 2*STRUCT_FACTORS_PRODS(Q,I,J)*SIN(ARG)/ARG
               END DO
            END DO
         END DO

         !apply logarithms
         LOGI(1:NQPOINTS) = LOG10(I1(1:NQPOINTS))
         !Define target curve
         I0(1:NQPOINTS) = 10**TARGET_CURVE(1:NQPOINTS)
         !Use correction
         IF (USE_MEAN_CORRECTION) THEN
            LOGI(1:NQPOINTS) = LOGI(1:NQPOINTS) + MEAN_CORR_CURVE(1:NQPOINTS)
            I1(1:NQPOINTS) = 10**LOGI(1:NQPOINTS)
         END IF

         ! Compute SAXS E and F, we try to use memorisation as much as possible
         CSCALE = I0(1)/I1(1)
         FSCALE = 1.0D0/(NQPOINTS*I1(1)**2)
         EDUMMY = 0.0D0
         DO Q=2,NQPOINTS
            QPHYS = QPOINTS(Q)
            DIFF = CSCALE*I1(Q) - I0(Q)
            DIFFQ = DIFF*QPHYS
            EDUMMY = EDUMMY + DIFFQ**2
            F_PRE = -2*QPHYS*DIFFQ*FSCALE
            DO I = 1,NPARTICLES
               IDX = 3*I
               DO J = I+1,NUM_TOT
                  JDX = 3*J
                  R(1:3) = VEC_DIFF(X_TOT(IDX-2:IDX),X_TOT(JDX-2:JDX))
                  R2 = DOT_PRODUCT(R,R)
                  QR = QPHYS*SQRT(R2)
                  FGRAIN = (R/R2)*STRUCT_FACTORS_PRODS(Q,I,J)*(COS(QR)-SIN(QR)/QR)
                  FACTING = 2*FGRAIN*F_PRE
                  FSAXS(IDX-2:IDX) = FSAXS(IDX-2:IDX) + FACTING
                  IF (J.LT.NPARTICLES) THEN
                     FSAXS(JDX-2:JDX) = FSAXS(JDX-2:JDX) - FACTING
                  END IF
               END DO
            END DO
         END DO
         !normalise
         ESAXS = EDUMMY/(NQPOINTS*I1(1))

         ! deallocate arrays
         DEALLOCATE(X_TOT,ALLNAMES,STRUCT_FACTORS,DISTMATRIX,STRUCT_FACTORS_PRODS)
         IF (ALLOCATED(XSOL)) DEALLOCATE(XSOL)
      END SUBROUTINE GENERATE_SAXS_CURVE_HYDRATE

      SUBROUTINE INIT_HASH_CG()
         USE HASH_INT
         IMPLICIT NONE
         INTEGER :: FUNIT, STATUS, GRAINIDX, GRAIN_HASH_NUM, Q
         CHARACTER(LEN=150) :: PARAMETER_LINE, DATFILE
         CHARACTER(LEN=5) :: BUFFERKEY
         REAL(KIND=REAL64) :: BUFFER1, BUFFER2, BUFFER3, BUFFER4, TEMP

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
         !REMOVE after benchmarking (needed to run both sets of calculations)
         !CALL HASH_DESTROY()
         CALL HASH_INIT()

         ! second parse - get grain names hashed
         GRAINIDX=0
         FUNIT = GETUNIT()
         OPEN(UNIT=FUNIT, FILE=PARAMETER_FILE, STATUS='OLD')
         DO WHILE (GRAINIDX.LT.NGRAINS)
            READ(FUNIT, '(A)', IOSTAT=STATUS) PARAMETER_LINE
            IF (STATUS.NE.0) EXIT
            IF (PARAMETER_LINE(1:1).NE.COMMENT_PATTERN) THEN
               GRAINIDX = GRAINIDX + 1
               READ(PARAMETER_LINE,'(A5)') BUFFERKEY
               CALL HASH_SET(BUFFERKEY,GRAINIDX)
               GRAINS(GRAINIDX) = BUFFERKEY
            END IF
         END DO
         CLOSE(FUNIT)

         ! Fill arrays with grain information
         GRAINIDX=0
         FUNIT = GETUNIT()
         OPEN(UNIT=FUNIT, FILE=PARAMETER_FILE, STATUS='OLD')
         DO WHILE (GRAINIDX.LT.NGRAINS)
            READ(FUNIT, '(A)', IOSTAT=STATUS) PARAMETER_LINE
            IF (STATUS.NE.0) EXIT
            IF (PARAMETER_LINE(1:1).NE.COMMENT_PATTERN) THEN
               GRAINIDX = GRAINIDX + 1
               READ(PARAMETER_LINE,'(A5,2F4.1,2F6.1,1X,A1)') &
                    BUFFERKEY, BUFFER1, BUFFER2, BUFFER3, BUFFER4
               CALL HASH_GET(BUFFERKEY,GRAIN_HASH_NUM)
               GRAINS_EXCLR(GRAIN_HASH_NUM) = BUFFER1
               GRAINS_CUTOFF(GRAIN_HASH_NUM) = BUFFER2
               GRAINS_EXCLV(GRAIN_HASH_NUM) = BUFFER3
               GRAINS_EDENSITY(GRAIN_HASH_NUM) = BUFFER4
               IF (BUFFERKEY.EQ.(SOLVENT_NAME(1:3)//"o")) THEN
                  CALL HASH_GET(SOLVENT_NAME,GRAIN_HASH_NUM)
                  IF (GRAIN_HASH_NUM.EQ.0) THEN
                     WRITE(*,*) " init_hash_cg> Must put solvent grain before outer solvent shell"
                     STOP
                  END IF
                  SOLVENT_CONTRAST = BUFFER4/GRAINS_EDENSITY(GRAIN_HASH_NUM)
               END IF
            END IF
         END DO
         CLOSE(FUNIT)

         ! Build the form factors
         DO GRAINIDX=1,NGRAINS
            IF (IN_SOLUTION_CURVE.AND.(.NOT.EXPLICIT_SOL_CONTR)) THEN
               DATFILE = 'dat/ff' // TRIM(GRAINS(GRAINIDX)) // '.cor.awk.dat'
            ELSE
               DATFILE = 'dat/ff' // TRIM(GRAINS(GRAINIDX)) // '.awk.dat'
            END IF
            FUNIT = GETUNIT()
            OPEN(UNIT=FUNIT, FILE=DATFILE, STATUS='OLD')
            DO Q=1,NQPOINTS
               READ(FUNIT, '(F8.3,F10.8)') TEMP, FACTOR_GRAINS(Q,GRAINIDX)
            END DO
            CLOSE(FUNIT)
            IF (TRIM(GRAINS(GRAINIDX)).EQ.TRIM(SOLVENT_NAME)) THEN
               FACTOR_GRAINS(:,GRAINIDX) = SAXS_W_SHELL * FACTOR_GRAINS(:,GRAINIDX)
            END IF
         END DO

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

      SUBROUTINE FILL_STRUCTURE_FACTOR_CG(NPARTS, PARTNAMES, STRUCT_F, STRUCT_F_PRODS)
         USE HASH_INT
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NPARTS ! number of particles (either NPARTICLES or NUM_TOT)
         CHARACTER(LEN=4), INTENT(IN) :: PARTNAMES(NPARTS)
         REAL(KIND=REAL64), INTENT(OUT) :: STRUCT_F(NQPOINTS,NPARTS)
         REAL(KIND=REAL64), INTENT(OUT) :: STRUCT_F_PRODS(NQPOINTS,NPARTS,NPARTS)
         INTEGER :: GRAIN_HASH_NUM, IDX, I, J, K

         STRUCT_F(:,:) = 0.0
         DO IDX=1,NPARTS
            CALL HASH_GET(ADJUSTL(TRIM(PARTNAMES(IDX))), GRAIN_HASH_NUM)
            IF (GRAIN_HASH_NUM.NE.0) THEN
               STRUCT_F(1:NQPOINTS,IDX) = FACTOR_GRAINS(1:NQPOINTS,GRAIN_HASH_NUM)
            ELSE
               WRITE(*,'(A)') " fill_struct_factor_cg> Could not find structure factor parameters for Grain: "
               WRITE(*,'(I6,A)') IDX, PARTNAMES(IDX)
               STOP 
            END IF        
         END DO

         ! get products of structure factors
         STRUCT_F_PRODS(:,:,:) = 0.0D0
         DO I=1,NPARTS
            DO J = I,NPARTS
               DO K=1,NQPOINTS
                  STRUCT_F_PRODS(K,I,J) = STRUCT_F(K,I)*STRUCT_F(K,J)
                  STRUCT_F_PRODS(K,J,I) = STRUCT_F_PRODS(K,I,J)
               END DO
            END DO
         END DO
      END SUBROUTINE FILL_STRUCTURE_FACTOR_CG

      SUBROUTINE FILL_HALF_DIST_MATRIX(NUM_GRAINS,X,DMATRIX)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NUM_GRAINS
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NUM_GRAINS)
         REAL(KIND=REAL64), INTENT(OUT) :: DMATRIX(NUM_GRAINS,NUM_GRAINS)
         INTEGER :: I, J, IDX, JDX
         REAL(KIND=REAL64) :: RIJ

         DMATRIX(:,:) = 0.0D0
         DO I = 1,NUM_GRAINS-1
            IDX = 3*I
            DO J = I+1,NUM_GRAINS
               JDX = 3*J
               RIJ = SQRT((X(IDX-2)-X(JDX-2))**2 + (X(IDX-1)-X(JDX-1))**2 +(X(IDX)-X(JDX))**2)
               DMATRIX(I,J) = RIJ
               DMATRIX(J,I) = RIJ
            END DO
         END DO

      END SUBROUTINE FILL_HALF_DIST_MATRIX

      SUBROUTINE HYDRATE_STRUCTURE(XYZ)
         USE VAR_DEFS, ONLY: ATNAMES => IGRAPH
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: XYZ(NPARTICLES)
         REAL(KIND=REAL64) :: CUTOFF_INTERNAL(NPARTICLES), CUTOFF_EXTERNAL(NPARTICLES)
         REAL(KIND=REAL64) :: MINMAXBOX(6)
         INTEGER :: XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN, IX, IY, IZ  !variables for iteration over box dimensions
         INTEGER :: NMAXSOL = 5000 ! vairbale used to assign local arrays
         REAL(KIND=REAL64), ALLOCATABLE :: XSOL_DUMMY(:), DUMMY(:)
         INTEGER :: NSOL ! actual number of solvent grains
         INTEGER :: CONTACT_NUM !number of contacts with grains
         REAL(KIND=REAL64) :: X, Y, Z !grid point position
         INTEGER :: I, IDX 
         REAL(KIND=REAL64) :: CUTINT, CUTEXT, RIJ !grain specific cutoff and distances

         ! get cutoffs for all grains
         CALL FILL_CUTOFFS(ATNAMES,CUTOFF_INTERNAL,CUTOFF_EXTERNAL)
         ! get box dimensions
         CALL BOX_DIMENSIONS(MINMAXBOX,SQRT(MAXVAL(CUTOFF_EXTERNAL)))
         ! fromt his get the number of points for x, y and z for iteration
         XMAX = INT(MINMAXBOX(1)/DX)
         XMIN = INT(MINMAXBOX(2)/DX)
         YMAX = INT(MINMAXBOX(3)/DX)
         YMIN = INT(MINMAXBOX(4)/DX)
         ZMAX = INT(MINMAXBOX(5)/DX)
         ZMIN = INT(MINMAXBOX(6)/DX)

         NSOL = 0
         ALLOCATE(XSOL_DUMMY(NMAXSOL))
         XSOL_DUMMY=0.0D0
         ! iterate over all grid points and place grains
         DO IX=XMIN,XMAX
            X = IX*DX
            DO IY=YMIN,YMAX
               Y = IY*DX
               DO IZ=ZMIN,ZMAX
                  Z = IZ*DX
                  CONTACT_NUM = 0 
                  ! cycle over all grains - are there clashes or good contacts?
                  DO I=1,NPARTICLES
                     CUTINT = CUTOFF_INTERNAL(I)
                     CUTEXT = CUTOFF_EXTERNAL(I)
                     IDX = 3*I
                     !get distance between grain and grid point
                     RIJ = (XYZ(IDX-2)-X)**2 
                     IF (RIJ.LE.CUTEXT) THEN
                        RIJ = RIJ + (XYZ(IDX-1)-Y)**2
                        IF (RIJ.LE.CUTEXT) THEN
                           RIJ = RIJ + (XYZ(IDX)-Z)**2
                           ! detect if there is a clash
                           IF (RIJ.LE.CUTINT) THEN
                              CONTACT_NUM = 0
                              EXIT
                           ! if not, check whether we are within the outer cutoff
                           ELSE IF (RIJ.LE.CUTEXT) THEN
                              CONTACT_NUM = CONTACT_NUM + 1
                           END IF
                        END IF
                     END IF
                  END DO
                  IF (CONTACT_NUM.GT.0) THEN
                     NSOL = NSOL+1
                     IF (NSOL.GT.NMAXSOL) THEN
                        ALLOCATE(DUMMY(NMAXSOL))
                        DUMMY(:) = XSOL_DUMMY(:)
                        DEALLOCATE(XSOL_DUMMY)
                        ALLOCATE(XSOL_DUMMY(2*NMAXSOL))
                        XSOL_DUMMY(1:NMAXSOL) = DUMMY(NMAXSOL)
                        NMAXSOL = 2*NMAXSOL
                        DEALLOCATE(DUMMY)
                     END IF
                     XSOL_DUMMY(3*NSOL-2) = X
                     XSOL_DUMMY(3*NSOL-1) = X
                     XSOL_DUMMY(3*NSOL) = X
                  END IF
               END DO
            END DO
         END DO
         NUM_SOL = NSOL
         ALLOCATE(XSOL(NUM_SOL))
         XSOL(1:NUM_SOL) = XSOL_DUMMY(1:NUM_SOL)
         DEALLOCATE(XSOL_DUMMY)
      END SUBROUTINE HYDRATE_STRUCTURE

      SUBROUTINE FILL_CUTOFFS(PARTNAMES, CUTOFF_INTERNAL, CUTOFF_EXTERNAL)
         USE HASH_INT
         IMPLICIT NONE
         CHARACTER(LEN=4), INTENT(IN) :: PARTNAMES(NPARTICLES)
         REAL(KIND=REAL64), INTENT(OUT) :: CUTOFF_INTERNAL(NPARTICLES)
         REAL(KIND=REAL64), INTENT(OUT) :: CUTOFF_EXTERNAL(NPARTICLES)   
         INTEGER :: GRAIN_HASH_NUM, IDX  
         
         DO IDX=1,NPARTICLES
            CALL HASH_GET(PARTNAMES(IDX),GRAIN_HASH_NUM)
            IF (GRAIN_HASH_NUM.EQ.0) THEN
               WRITE(*,'(A)') " fill_struct_factor_cg> Could not find structure factor parameters for Grain: "
               WRITE(*,'(I6,A)') IDX, PARTNAMES(IDX)
               STOP 
            END IF
            CUTOFF_INTERNAL(IDX) = (0.50 * DX * GRAINS_CUTOFF(GRAIN_HASH_NUM))**2
            CUTOFF_EXTERNAL(IDX) = (0.99 * (N_SHELLS + 0.5) * DX * GRAINS_CUTOFF(GRAIN_HASH_NUM))**2
         END DO

      END SUBROUTINE FILL_CUTOFFS


      SUBROUTINE BOX_DIMENSIONS(MINMAXBOX, DPADDING)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(OUT) :: MINMAXBOX(6)
         REAL(KIND=REAL64), INTENT(IN) :: DPADDING
         ! PDB is limited to 3 digits before the dot, so [-1000..1000] is the range
         MINMAXBOX=(/1000.0,-1000.0,1000.0,-1000.0,1000.0,-1000.0/)
         ! Pad by padding_length nm in every direction, just to be safe
         ! Warning : this value should be modified for really big coarse grained solvent molecule
         MINMAXBOX(1) = MINMAXBOX(1) + DPADDING
         MINMAXBOX(2) = MINMAXBOX(2) - DPADDING
         MINMAXBOX(3) = MINMAXBOX(3) + DPADDING
         MINMAXBOX(4) = MINMAXBOX(4) - DPADDING
         MINMAXBOX(5) = MINMAXBOX(5) + DPADDING
         MINMAXBOX(6) = MINMAXBOX(6) - DPADDING         
      END SUBROUTINE BOX_DIMENSIONS

END MODULE SAXS_CALCS
