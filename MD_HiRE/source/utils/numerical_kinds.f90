!> @file
!> Contains NUMKIND module with definitions of numerical kinds

!> Module to provide consistent kinds for floats and integers
!> @brief
MODULE NUMKIND
   !> integer 32-bit
   INTEGER, PARAMETER  :: INT32  = SELECTED_INT_KIND(9)
   !> integer 64-bit
   INTEGER, PARAMETER  :: INT64  = SELECTED_INT_KIND(18)
   !> real, double precision 32-bit
   INTEGER, PARAMETER  :: REAL32 = SELECTED_REAL_KIND(6, 37)
   !> real, double precision 64-bit
   INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)

END MODULE NUMKIND
