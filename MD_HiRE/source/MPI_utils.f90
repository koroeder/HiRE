MODULE MD_MPI
   USE NUMKIND
   IMPLICIT NONE

   CONTAINS

      SUBROUTINE COMMUNICATE_SETTINGS()
         USE MD_COMMONS

         ! use BCAST for globals

         ! use SEND/RECV for temperature and lambda
         
      END SUBROUTINE COMMUNICATE_SETTINGS


END MODULE MD_MPI