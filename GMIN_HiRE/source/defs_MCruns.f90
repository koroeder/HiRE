MODULE DEFS_MCRUNS
! Module containing the definitions used in the BH routines and minimisers that used to be COMMON blocks
  USE PREC
  IMPLICIT NONE

  !Test flags
  LOGICAL :: FTEST
  LOGICAL :: STUCK
  INTEGER :: NQTOT, ITER, POINT
  INTEGER :: ISPT, IYPT, NPT
  REAL(KIND = REAL64) :: QSTART, QFINISH
  REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: DIAG, W
END MODULE DEFS_MCRUNS
