MODULE MOD_VEL
   USE NUMKIND
   USE MD_COMMONS, ONLY: NATOMS, MASSES, MYUNIT
   IMPLICIT NONE

   CONTAINS
      
      SUBROUTINE CREATE_RND_VEL(TEMP, COORDS, VEL, COM, P_COM)
         USE RAND_ROUTINES, ONLY: RAND_NORMAL
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: TEMP             ! temperature
         REAL(KIND=REAL64), INTENT(IN) :: COORDS(3*NATOMS) ! coordinates
         REAL(KIND=REAL64), INTENT(OUT) :: VEL(3*NATOMS)   ! velocity
         REAL(KIND=REAL64), INTENT(OUT) :: COM(3)          ! centre of mass
         REAL(KIND=REAL64), INTENT(OUT) :: P_COM(3)        ! momentum of centre of mass
         INTEGER :: I, J
         REAL(KIND=REAL64), PARAMETER :: STDEV = 1.0D0
         REAL(KIND=REAL64), PARAMETER :: MEAN = 0.0D0
         REAL(KIND=REAL64) :: ATMASS, CURRVEL, TOTALMASS
         
         TOTALMASS = SUM(MASSES)
         COM(1:3) = 0.0D0
         P_COM(1:3) = 0.0D0
         DO I=1,NATOMS
            ATMASS = MASSES(I)
            DO J=1,3
               CALL RAND_NORMAL(STDEV,MEAN,CURRVEL)
               VEL(3*(I-1)+J) = CURRVEL
               COM(J) = COM(J) + COORDS(3*(I-1)+J)*ATMASS
               P_COM(J) = P_COM(J) + CURRVEL*ATMASS
            ENDDO
         END DO
      END SUBROUTINE CREATE_RND_VEL


      SUBROUTINE THERMALISE()

         !from thermalize routine

         ! 1. remove angular momentum and rotations
         ! 2. cycle over thermalisation steps
         !   a) rescale velocity
         !   b) cycle of equilibration steps
         !      i) potentially remove total lin mom and ang mom
         !      ii) potentially rescale velocities
         !      iii) velocity Verlet step
         !   c) report statistics for this

         ! rescaling velocity 2 (a)
         ! 1. get current temperature (kinetic energy divided by #dof)
         ! 2. rescale with scale = sqrt(target temp/current temp)
         !    vel = vel*scale
      END SUBROUTINE THERMALISE
END MODULE MOD_VEL