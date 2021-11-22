MODULE NAparams
  USE PREC_HIRE, ONLY: REAL64
  IMPLICIT NONE

  INTEGER, PARAMETER :: SCORESIZE=47
  REAL(KIND = REAL64) :: SCORE_RNA(SCORESIZE)
  
  !RNA bases
  INTEGER, ALLOCATABLE :: BLIST(:)
  INTEGER, ALLOCATABLE :: BTYPE(:)
  INTEGER, ALLOCATABLE :: BPROT(:)   
  INTEGER, ALLOCATABLE :: BOCC(:)
  REAL(KIND = REAL64), ALLOCATABLE :: BPCH(:)
  LOGICAL, ALLOCATABLE :: BP_CURR(:,:)
  
  !Cut off definitions
  REAL(KIND = REAL64), PARAMETER :: CUT   = 1.0D2 ! NO CUTOFF
  REAL(KIND = REAL64), PARAMETER :: SCNB  = 8.0D1 ! divide the 1-4 VDW interactions by 8.0
  REAL(KIND = REAL64), PARAMETER :: SCEE  = 8.0D1 ! divide the 1-4 ELEC interactions by 2.0
  REAL(KIND = REAL64), PARAMETER :: IDIEL = 0.0D0 ! 2.0 dielectric constant epsilon = 2r
  REAL(KIND = REAL64), PARAMETER :: DIELC = 1.0D0 ! Dielectric constant
           
  REAL(KIND = REAL64), PARAMETER :: Z = 0.0D0
  
  CHARACTER(LEN=1), DIMENSION(4) :: basenum = (/'G','A','C','U'/)        !
  
  INTEGER :: NPARAM(4,4) !
  
  REAL(KIND = REAL64) :: wc                 !rna/dna parameter: wc coef
  REAL(KIND = REAL64) :: wcCanonic          !rna/dna parameter: wcCanonic coef  
  REAL(KIND = REAL64) :: noWc               !rna/dna parameter: noWc coef
  REAL(KIND = REAL64) :: tit                !rna parameter titration
  REAL(KIND = REAL64) :: noWCq              !rna parameter hb charged base
   
  !Variables for titration
  INTEGER, ALLOCATABLE :: fatortit(:)       !titration start
  
  !RNA cutoffs
  REAL(KIND = REAL64), PARAMETER :: rcut2_caca_scsc_out = 1.44D2
  REAL(KIND = REAL64), PARAMETER :: rcut2_caca_scsc_in = 1.00D2
  REAL(KIND = REAL64), PARAMETER :: rcut2_hb_mcmc_out = 1.44D2
  REAL(KIND = REAL64), PARAMETER :: rcut2_hb_mcmc_in = 1.00D2
  REAL(KIND = REAL64), PARAMETER :: rcut2_4b_out = 1.44D2
  REAL(KIND = REAL64), PARAMETER :: rcut2_4b_in = 1.00D2
  REAL(KIND = REAL64), PARAMETER :: rcut2_lj_out = 1.44D2
  REAL(KIND = REAL64), PARAMETER :: rcut2_lj_in = 1.00D2
  
  SAVE
END MODULE NAparams

MODULE RNA_HB_PARAMS
  USE PREC_HIRE, ONLY: REAL64
  IMPLICIT NONE
  INTEGER :: NPARAM(4,4) !
  REAL(KIND = REAL64) :: alpa(4,4)            ! 
  REAL(KIND = REAL64) :: alpb(4,4)            !

  REAL(KIND = REAL64) :: alpam(6,4,4)         ! 
  REAL(KIND = REAL64) :: alpbm(6,4,4)         !
  REAL(KIND = REAL64) :: calpam(6,4,4)        ! 
  REAL(KIND = REAL64) :: calpbm(6,4,4)        !
  REAL(KIND = REAL64) :: salpam(6,4,4)        ! 
  REAL(KIND = REAL64) :: salpbm(6,4,4)        !
  REAL(KIND = REAL64) :: dREF(6,4,4)          ! 
 
  REAL(KIND = REAL64) :: s(6,4,4,2,2)         !
  
  REAL(KIND = REAL64) :: planarityDistEq(4,3) !
  
  CONTAINS
   
   SUBROUTINE FILL_RNA_HB_PARAMS
      USE NAparams, ONLY: WC, WCCanonic, noWC, TIT, noWCq, Z   
      USE NUM_DEFS, ONLY: PI
      
      alpam(:,:,:) = 2.0D0
      alpbm(:,:,:) = 2.0D0
      !         A-A
      dREF(1:3,2,2) =  (/ 5.63, 6.84, 5.92/)       !WWt, HHt, HWt                    
      alpam(1:3,2,2) = (/ 2.40, 1.05, 2.72/)
      alpbm(1:3,2,2) = (/ 2.40, 1.08, 1.29/)
      s(1:3,2,2,1,1) = (/2*wc, 2*noWc, 2*noWc/)    !Br2  1,1 -> q1=0, q2=0
      s(1:3,2,2,1,2) = (/    z,2*noWcq, z/)        !Br2  1,2 -> q1=0, q2=1
      s(1:3,2,2,2,1) = (/    z,2*noWcq, 2*noWcq/)  !Br2  2,1 -> q1=1, q2=0
      s(1:3,2,2,2,2) = (/    z,2*noWcq, z/)        !Br2  2,2 -> q1=1, q2=1
      Nparam(2,2) = 3  !Br3  CHECK angle sign when order is inverted when one is negative

      !         A-C
      dREF(1:3,2,3) = (/7.26, 5.78, 4.78/)         !wsc, HWt, +Wc  
      alpam(1:3,2,3) = (/2.36, 1.36, 2.64/)
      alpbm(1:3,2,3) = (/0.75, 2.30, 1.78/)   
      s(1:3,2,3,1,1) = (/ 2*noWc,  2*noWc, z /)    !Br2  1,1 -> q1=0, q2=0
      s(1:3,2,3,1,2) = (/ 2*noWcq, z, z/)          !Br2  1,2 -> q1=0, q2=1
      s(1:3,2,3,2,1) = (/ z, 2*noWcq, 2*tit/)      !Br2  2,1 -> q1=1, q2=0
      s(1:3,2,3,2,2) = (/ z, z, z/)                !Br2  2,2 -> q1=1, q2=1

      Nparam(2,3) = 3
      dREF(:,3,2) = dREF(:,2,3)
      alpam(:,3,2) = alpbm(:,2,3)
      alpbm(:,3,2) = alpam(:,2,3)
      s(:,3,2,1,1) = s(:,2,3,1,1)
      s(:,3,2,1,2) = s(:,2,3,2,1)
      s(:,3,2,2,1) = s(:,2,3,1,2)
      s(:,3,2,2,2) = s(:,2,3,2,2)
      Nparam(3,2) = Nparam(2,3)

      !         A-G
      dREF(1:4,2,1) = (/ 4.88, 6.63, 6.17, 6.05/)         !  WW_c, HSt, sst, +Hc  
      alpam(1:4,2,1) = (/ 3.04, 1.19, -2.01, 2.62/)
      alpbm(1:4,2,1) = (/ 2.58, -1.69, -1.57, 0.89/)
      s(1:4,2,1,1,1) = (/ 2*wc, 2*noWc, 2*noWc, z/)       !Br2  1,1 -> q1=0, q2=0
      s(1:4,2,1,1,2) = (/ z, 2*noWcq, 2*noWcq, z/)        !Br2  1,2 -> q1=0, q2=1
      s(1:4,2,1,2,1) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)    !Br2  2,1 -> q1=1, q2=0
      s(1:4,2,1,2,2) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)    !Br2  2,2 -> q1=1, q2=1

      Nparam(2,1) = 4
      dREF(:,1,2) = dREF(:,2,1)
      alpam(:,1,2) = alpbm(:,2,1)
      alpbm(:,1,2) = alpam(:,2,1)
      s(:,1,2,1,1) = s(:,2,1,1,1)
      s(:,1,2,1,2) = s(:,2,1,2,1)
      s(:,1,2,2,1) = s(:,2,1,1,2)
      s(:,1,2,2,2) = s(:,2,1,2,2)
      Nparam(1,2) = Nparam(2,1)

      !         A-U
      dREF(1:3,2,4) = (/ 4.92, 6.0, 6.0/)        ! WWc, HWt, HWc  5.78, 5.89
      alpam(1:3,2,4) = (/ 2.84, 0.91, 0.93/)
      alpbm(1:3,2,4) = (/ 2.36, 1.82, 2.46/)
      s(1:3,2,4,1,1) = (/ 2.2*wcCanonic, 2*noWc, 2*noWc/) !Br2  1,1 -> q1=0, q2=0 
                                          !OKKIO : cambiato a mano 14.4 -> 16!!!
      s(1:3,2,4,1,2) = (/ z, z, z/)                       !Br2  1,2 -> q1=0, q2=1
      s(1:3,2,4,2,1) = (/ z, 2*noWcq, 2*noWcq/)           !Br2  2,1 -> q1=1, q2=0
      s(1:3,2,4,2,2) = (/ z, z, z/)                       !Br2  2,2 -> q1=1, q2=1
      Nparam(2,4) = 3
          
      dREF(:,4,2) = dREF(:,2,4)
      alpam(:,4,2) = alpbm(:,2,4)
      alpbm(:,4,2) = alpam(:,2,4)
      s(:,4,2,1,1) = s(:,2,4,1,1)
      s(:,4,2,1,2) = s(:,2,4,2,1)
      s(:,4,2,2,1) = s(:,2,4,1,2)
      s(:,4,2,2,2) = s(:,2,4,2,2)
      Nparam(4,2) = Nparam(2,4)

! !         C-C
!           dREF(1:1,3,3) =  (/ 4.91/)  !Br3 NO significant CC pairing
!           alpam(1:1,3,3) = (/ 2.22/)
!           alpbm(1:1,3,3) = (/ 2.24/)
!           s(1:1,3,3,1,1) = (/ z /) !Br2  1,1 -> q1=0, q2=0
!           s(1:1,3,3,1,2) = (/ z /) !Br2  1,2 -> q1=0, q2=1
!           s(1:1,3,3,2,1) = (/ z /) !Br2  2,1 -> q1=1, q2=0
!           s(1:1,3,3,2,2) = (/ z /) !Br2  2,2 -> q1=1, q2=1
!           Nparam(3,3) = 1

      !         C-G
      dREF(1:3,3,1) = (/ 4.75, 5.28, 5.68/)              ! WWc, WWt, +Wc
      alpam(1:3,3,1) = (/ 2.17, 1.64, 1.92/)
      alpbm(1:3,3,1) = (/ 2.71, -3.07, 2.34/)
      s(1:3,3,1,1,1) = (/ 2.6*wcCanonic, 2*wc, z/)       !Br2  1,1 -> q1=0, q2=0
      s(1:3,3,1,1,2) = (/ z, z, z/)                      !Br2  1,2 -> q1=0, q2=1
      s(1:3,3,1,2,1) = (/ z, z, 1*tit/)                  !Br2  2,1 -> q1=1, q2=0
      s(1:3,3,1,2,2) = (/ z, z, z/)                      !Br2  2,2 -> q1=1, q2=1
      Nparam(3,1) = 3
          
      dREF(:,1,3) = dREF(:,3,1)
      alpam(:,1,3) = alpbm(:,3,1)
      alpbm(:,1,3) = alpam(:,3,1)
      s(:,1,3,1,1) = s(:,3,1,1,1)
      s(:,1,3,1,2) = s(:,3,1,2,1)
      s(:,1,3,2,1) = s(:,3,1,1,2)
      s(:,1,3,2,2) = s(:,3,1,2,2)
      Nparam(1,3) = Nparam(3,1)

      !         C-U
      dREF(1:1,3,4) = (/ 4.81 /)             ! WWc
      alpam(1:1,3,4) = (/ 2.02 /)
      alpbm(1:1,3,4) = (/-2.50 /)            !Br3 CHECK when particles are inverted
      s(1:1,3,4,1,1) = (/2*wc/)              !Br2  1,1 -> q1=0, q2=0     
      s(1:1,3,4,1,2) = (/ z /)               !Br2  1,2 -> q1=0, q2=1
      s(1:1,3,4,2,1) = (/ z /)               !Br2  2,1 -> q1=1, q2=0
      s(1:1,3,4,2,2) = (/ z /)               !Br2  2,2 -> q1=1, q2=1
      Nparam(3,4) = 1
          
      dREF(:,4,3) = dREF(:,3,4)
      alpam(:,4,3) = alpbm(:,3,4)
      alpbm(:,4,3) = alpam(:,3,4)
      s(:,4,3,1,1) = s(:,3,4,1,1)
      s(:,4,3,1,2) = s(:,3,4,2,1)
      s(:,4,3,2,1) = s(:,3,4,1,2)
      s(:,4,3,2,2) = s(:,3,4,2,2)
      Nparam(4,3) = Nparam(3,4)

      !         G-G
      dREF(1:3,1,1) =  (/ 6.00, 6.00, 6.73 /)        !     HWc, HWt, SSt   6.22 -> 6.0, 6.25 -> 6.00
      alpam(1:3,1,1) = (/ 1.27, 3.02, -1.78 /)
      alpbm(1:3,1,1) = (/ 2.90, 0.82, -1.83 /)
      s(1:3,1,1,1,1) = (/ 2*noWc, 2*noWc, 2*noWc /)      !Br2  1,1 -> q1=0, q2=0     
      s(1:3,1,1,1,2) = (/ z, z, 2*noWcq /)               !Br2  1,2 -> q1=0, q2=1
      s(1:3,1,1,2,1) = (/  2*noWcq, 2*noWcq, 2*noWcq/)   !Br2  2,1 -> q1=1, q2=0
      s(1:3,1,1,2,2) = (/ z, z, 2*noWcq/)                !Br2  2,2 -> q1=1, q2=1
      Nparam(1,1) = 3

      !         G-U
      dREF(1:1,1,4) = (/ 5.05 /)               ! WWc         
      alpam(1:1,1,4) = (/ 2.29 /)
      alpbm(1:1,1,4) = (/ 1.68 /)
      s(1:1,1,4,1,1) = (/ 2.1*wcCanonic /)     !Br2  1,1 -> q1=0, q2=0 !OKKIO : cambiato a mano 14.7 -> 16!!!
      s(1:1,1,4,1,2) = (/ z /)                 !Br2  1,2 -> q1=0, q2=1
      s(1:1,1,4,2,1) = (/ z /)                 !Br2  2,1 -> q1=1, q2=0
      s(1:1,1,4,2,2) = (/ z /)                 !Br2  2,2 -> q1=1, q2=1
      Nparam(1,4) = 1
          
      dREF(:,4,1) = dREF(:,1,4)
      alpam(:,4,1) = alpbm(:,1,4)
      alpbm(:,4,1) = alpam(:,1,4)
      s(:,4,1,1,1) = s(:,1,4,1,1)
      s(:,4,1,1,2) = s(:,1,4,2,1)
      s(:,4,1,2,1) = s(:,1,4,1,2)
      s(:,4,1,2,2) = s(:,1,4,2,2)
      Nparam(4,1) = Nparam(1,4)

      !         U-U
      dREF(1:3,4,4) =  (/ 4.94, 4.84, 5.63 /)           ! WWc, WWt, wht
      alpam(1:3,4,4) = (/ 1.85, -1.88, 2.36 /)
      alpbm(1:3,4,4) = (/ 2.48, 1.71, -2.57 /)
      s(1:3,4,4,1,1) = (/ 2*wc, 2*wc, 2*noWc /)        !Br2  1,1 -> q1=0, q2=0
      s(1:3,4,4,1,2) = (/ z, z, 2*noWcq /)                !Br2  1,2 -> q1=0, q2=1
      s(1:3,4,4,2,1) = (/ z, z, z /)                    !Br2  2,1 -> q1=1, q2=0
      s(1:3,4,4,2,2) = (/ z, z, z /)                    !Br2  2,2 -> q1=1, q2=1
      Nparam(4,4) = 3

 
      !------------------------------------------------------------------
      !         G: CY G1 G2
      planarityDistEq(1,1:3) = (/ 0, 0, 0/)
      !         A: CY A1 A2
      planarityDistEq(2,1:3) = (/ 0, 0, 0/)
      !         C: CA CY C1
      planarityDistEq(3,1:3) = (/ 0, 0, 0/)
      !         U: CA CY U1
      planarityDistEq(4,1:3) = (/ 0, 0, 0/)  
      !------------------------------------------------------------------
      alpam = alpam+pi
      alpbm = alpbm+pi
      calpam = cos(alpam)
      salpam = sin(alpam)
      calpbm = cos(alpbm)
      salpbm = sin(alpbm)  
   END SUBROUTINE FILL_RNA_HB_PARAMS
END MODULE RNA_HB_PARAMS

MODULE DNA_HB_PARAMS
  USE PREC_HIRE, ONLY: REAL64
  IMPLICIT NONE
  INTEGER :: NPARAM(4,4) !
  REAL(KIND = REAL64) :: alpa(4,4)            ! 
  REAL(KIND = REAL64) :: alpb(4,4)            !

  REAL(KIND = REAL64) :: alpam(6,4,4)         ! 
  REAL(KIND = REAL64) :: alpbm(6,4,4)         !
  REAL(KIND = REAL64) :: calpam(6,4,4)        ! 
  REAL(KIND = REAL64) :: calpbm(6,4,4)        !
  REAL(KIND = REAL64) :: salpam(6,4,4)        ! 
  REAL(KIND = REAL64) :: salpbm(6,4,4)        !
  REAL(KIND = REAL64) :: dREF(6,4,4)          ! 
 
  REAL(KIND = REAL64) :: s(6,4,4,2,2)         !
  
  REAL(KIND = REAL64) :: planarityDistEq(4,3) !
    
  CONTAINS
   
   SUBROUTINE FILL_DNA_HB_PARAMS
      USE NAparams, ONLY: WC, WCCanonic, noWC, TIT, noWCq, Z
      USE NUM_DEFS, ONLY: PI
      
      alpam(:,:,:) = 2.0D0
      alpbm(:,:,:) = 2.0D0
      !         A-A
          dREF(1:3,2,2) =  (/ 5.63, 6.84, 5.92/)       !WWt, HHt, HWt                    
          alpam(1:3,2,2) = (/ 2.40, 1.05, 2.72/)
          alpbm(1:3,2,2) = (/ 2.40, 1.08, 1.29/)
          s(1:3,2,2,1,1) = (/2*wc, 2*noWc, 2*noWc/)     !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,2,1,2) = (/    z,2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,2,2,1) = (/    z,2*noWcq, 2*noWcq/)       !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,2,2,2) = (/    z,2*noWcq, z/)            !Br2  2,2 -> q1=1, q2=1
          Nparam(2,2) = 3                !Br3  CHECK angle sign when order is inverted when one is negative

!         A-C
          dREF(1:3,2,3) = (/7.26, 5.78, 4.78/)         !wsc, HWt, +Wc  
          alpam(1:3,2,3) = (/2.36, 1.36, 2.64/)
          alpbm(1:3,2,3) = (/0.75, 2.30, 1.78/)   
          s(1:3,2,3,1,1) = (/ 2*noWc,  2*noWc, z /) !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,3,1,2) = (/ 2*noWcq, z, z/)                 !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,3,2,1) = (/ z, 2*noWcq, 2.2*tit/)            !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,3,2,2) = (/ z, z, z/)                      !Br2  2,2 -> q1=1, q2=1

          Nparam(2,3) = 3
          dREF(:,3,2) = dREF(:,2,3)
          alpam(:,3,2) = alpbm(:,2,3)
          alpbm(:,3,2) = alpam(:,2,3)
          s(:,3,2,1,1) = s(:,2,3,1,1)
          s(:,3,2,1,2) = s(:,2,3,2,1)
          s(:,3,2,2,1) = s(:,2,3,1,2)
          s(:,3,2,2,2) = s(:,2,3,2,2)
          Nparam(3,2) = Nparam(2,3)

!         A-G
!          dREF(1:4,2,1) = (/ 4.88, 6.63, 6.17, 6.05/)         !  WW_c, HSt, sst, +Hc  
!          alpam(1:4,2,1) = (/ 3.04, 1.19, -2.01, 2.62/)
!          alpbm(1:4,2,1) = (/ 2.58, -1.69, -1.57, 0.89/)
!          s(1:4,2,1,1,1) = (/ 2*wc, 2*noWc, 2*noWc, z/)        !Br2  1,1 -> q1=0, q2=0
!          s(1:4,2,1,1,2) = (/ z, 2*noWcq, 2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
!          s(1:4,2,1,2,1) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)      !Br2  2,1 -> q1=1, q2=0
!          s(1:4,2,1,2,2) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)      !Br2  2,2 -> q1=1, q2=1

!          Nparam(2,1) = 4
!          dREF(:,1,2) = dREF(:,2,1)
!          alpam(:,1,2) = alpbm(:,2,1)
!          alpbm(:,1,2) = alpam(:,2,1)
!          s(:,1,2,1,1) = s(:,2,1,1,1)
!          s(:,1,2,1,2) = s(:,2,1,2,1)
!          s(:,1,2,2,1) = s(:,2,1,1,2)
!          s(:,1,2,2,2) = s(:,2,1,2,2)
          
!          Nparam(1,2) = Nparam(2,1)


!         A-G
          dREF(1:3,2,1) = (/ 4.88, 6.63, 6.05/)         !  WW_c, HSt, +Hc
          alpam(1:3,2,1) = (/ 3.04, 1.19,  2.62/)
          alpbm(1:3,2,1) = (/ 2.58, -1.69,  0.89/)
          s(1:3,2,1,1,1) = (/ 2*wc, 2*noWc, z/)        !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,1,1,2) = (/ z, 2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,1,2,1) = (/ z, 2*noWcq, 1.6*tit/)     !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,1,2,2) = (/ z, 2*noWcq, 1.6*tit/)      !Br2  2,2 -> q1=1, q2=1

          Nparam(2,1) = 3
          dREF(:,1,2) = dREF(:,2,1)
          alpam(:,1,2) = alpbm(:,2,1)
          alpbm(:,1,2) = alpam(:,2,1)
          s(:,1,2,1,1) = s(:,2,1,1,1)
          s(:,1,2,1,2) = s(:,2,1,2,1)
          s(:,1,2,2,1) = s(:,2,1,1,2)
          s(:,1,2,2,2) = s(:,2,1,2,2)

          Nparam(1,2) = Nparam(2,1)

!         A-U
          dREF(1:3,2,4) = (/ 4.92, 5.78, 5.89/)        ! WWc, HWt, HWc
          alpam(1:3,2,4) = (/ 2.84, 0.91, 0.93/)
          alpbm(1:3,2,4) = (/ 2.36, 1.82, 2.46/)
          s(1:3,2,4,1,1) = (/ 2.2*wcCanonic, 2*noWc, 2*noWc/)     !Br2  1,1 -> q1=0, q2=0      !OKKIO : cambiato a mano 14.4 -> 16!!!
          s(1:3,2,4,1,2) = (/ z, z, z/)                     !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,4,2,1) = (/ z, 2*noWcq, 2*noWcq/)          !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,4,2,2) = (/ z, z, z/)                     !Br2  2,2 -> q1=1, q2=1
          Nparam(2,4) = 3
          
          dREF(:,4,2) = dREF(:,2,4)
          alpam(:,4,2) = alpbm(:,2,4)
          alpbm(:,4,2) = alpam(:,2,4)
          s(:,4,2,1,1) = s(:,2,4,1,1)
          s(:,4,2,1,2) = s(:,2,4,2,1)
          s(:,4,2,2,1) = s(:,2,4,1,2)
          s(:,4,2,2,2) = s(:,2,4,2,2)
          Nparam(4,2) = Nparam(2,4)

! !         C-C
!           dREF(1:1,3,3) =  (/ 4.91/)  !Br3 NO significant CC pairing
!           alpam(1:1,3,3) = (/ 2.22/)
!           alpbm(1:1,3,3) = (/ 2.24/)
!           s(1:1,3,3,1,1) = (/ z /) !Br2  1,1 -> q1=0, q2=0
!           s(1:1,3,3,1,2) = (/ z /) !Br2  1,2 -> q1=0, q2=1
!           s(1:1,3,3,2,1) = (/ z /) !Br2  2,1 -> q1=1, q2=0
!           s(1:1,3,3,2,2) = (/ z /) !Br2  2,2 -> q1=1, q2=1
!           Nparam(3,3) = 1

!         C-G
          dREF(1:3,3,1) = (/ 4.75, 5.28, 5.68/)              ! WWc, WWt, +Wc
          alpam(1:3,3,1) = (/ 2.17, 1.64, 1.92/)
          alpbm(1:3,3,1) = (/ 2.71, -3.07, 2.34/)
          s(1:3,3,1,1,1) = (/ 2.6*wcCanonic, 2*wc, z/)       !Br2  1,1 -> q1=0, q2=0
          s(1:3,3,1,1,2) = (/ z, z, z/)                      !Br2  1,2 -> q1=0, q2=1
          s(1:3,3,1,2,1) = (/ z, z, 0.5*tit/)                  !Br2  2,1 -> q1=1, q2=0
          s(1:3,3,1,2,2) = (/ z, z, z/)                      !Br2  2,2 -> q1=1, q2=1
          Nparam(3,1) = 3
          
          dREF(:,1,3) = dREF(:,3,1)
          alpam(:,1,3) = alpbm(:,3,1)
          alpbm(:,1,3) = alpam(:,3,1)
          s(:,1,3,1,1) = s(:,3,1,1,1)
          s(:,1,3,1,2) = s(:,3,1,2,1)
          s(:,1,3,2,1) = s(:,3,1,1,2)
          s(:,1,3,2,2) = s(:,3,1,2,2)
          Nparam(1,3) = Nparam(3,1)

!         C-U
          dREF(1:1,3,4) = (/ 4.81 /)             ! WWc
          alpam(1:1,3,4) = (/ 2.02 /)
          alpbm(1:1,3,4) = (/-2.50 /)            !Br3 CHECK when particles are inverted
          s(1:1,3,4,1,1) = (/2*wc/)               !Br2  1,1 -> q1=0, q2=0     
          s(1:1,3,4,1,2) = (/ z /)                 !Br2  1,2 -> q1=0, q2=1
          s(1:1,3,4,2,1) = (/ z /)                 !Br2  2,1 -> q1=1, q2=0
          s(1:1,3,4,2,2) = (/ z /)                 !Br2  2,2 -> q1=1, q2=1
          Nparam(3,4) = 1
          
          dREF(:,4,3) = dREF(:,3,4)
          alpam(:,4,3) = alpbm(:,3,4)
          alpbm(:,4,3) = alpam(:,3,4)
          s(:,4,3,1,1) = s(:,3,4,1,1)
          s(:,4,3,1,2) = s(:,3,4,2,1)
          s(:,4,3,2,1) = s(:,3,4,1,2)
          s(:,4,3,2,2) = s(:,3,4,2,2)
          Nparam(4,3) = Nparam(3,4)

!         G-G
          dREF(1:3,1,1) =  (/ 6.22, 6.25, 6.73 /)        !     HWc, HWt, SSt
          alpam(1:3,1,1) = (/ 1.27, 3.02, -1.78 /)
          alpbm(1:3,1,1) = (/ 2.90, 0.82, -1.83 /)
          s(1:3,1,1,1,1) = (/ 2*noWc, 2*noWc, 2*noWc /)   !Br2  1,1 -> q1=0, q2=0     
          s(1:3,1,1,1,2) = (/ z, z, 2*noWcq /)             !Br2  1,2 -> q1=0, q2=1
          s(1:3,1,1,2,1) = (/  2*noWcq, 2*noWcq, 2*noWcq/)   !Br2  2,1 -> q1=1, q2=0
          s(1:3,1,1,2,2) = (/ z, z, 2*noWcq/)             !Br2  2,2 -> q1=1, q2=1
          Nparam(1,1) = 3

!         G-U
          dREF(1:1,1,4) = (/ 5.05 /)               ! WWc         
          alpam(1:1,1,4) = (/ 2.29 /)
          alpbm(1:1,1,4) = (/ 1.68 /)
          s(1:1,1,4,1,1) = (/ 2.1*wc /)            !Br2  1,1 -> q1=0, q2=0 !OKKIO : cambiato a mano 14.7 -> 16!!!
          s(1:1,1,4,1,2) = (/ z /)                 !Br2  1,2 -> q1=0, q2=1
          s(1:1,1,4,2,1) = (/ z /)                 !Br2  2,1 -> q1=1, q2=0
          s(1:1,1,4,2,2) = (/ z /)                 !Br2  2,2 -> q1=1, q2=1
          Nparam(1,4) = 1
          
          dREF(:,4,1) = dREF(:,1,4)
          alpam(:,4,1) = alpbm(:,1,4)
          alpbm(:,4,1) = alpam(:,1,4)
          s(:,4,1,1,1) = s(:,1,4,1,1)
          s(:,4,1,1,2) = s(:,1,4,2,1)
          s(:,4,1,2,1) = s(:,1,4,1,2)
          s(:,4,1,2,2) = s(:,1,4,2,2)
          Nparam(4,1) = Nparam(1,4)

!         U-U
          dREF(1:3,4,4) =  (/ 4.94, 4.84, 5.63 /)           ! WWc, WWt, wht
          alpam(1:3,4,4) = (/ 1.85, -1.88, 2.36 /)
          alpbm(1:3,4,4) = (/ 2.48, 1.71, -2.57 /)
          s(1:3,4,4,1,1) = (/ 2*wc, 2*wc, 2*noWc /)        !Br2  1,1 -> q1=0, q2=0
          s(1:3,4,4,1,2) = (/ z, z, 2*noWcq /)                !Br2  1,2 -> q1=0, q2=1
          s(1:3,4,4,2,1) = (/ z, z, z /)                    !Br2  2,1 -> q1=1, q2=0
          s(1:3,4,4,2,2) = (/ z, z, z /)                    !Br2  2,2 -> q1=1, q2=1
          Nparam(4,4) = 3


!     ------------------------------------------------------------------
!         G: CY G1 G2
!          planarityDistEq(1,1:3) = (/ 0, 0, 0/)
!         A: CY A1 A2
!          planarityDistEq(2,1:3) = (/ 0, 0, 0/)
!         C: CA CY C1
!          planarityDistEq(3,1:3) = (/ 0, 0, 0/)
!!         U: CA CY U1
!          planarityDistEq(4,1:3) = (/ 0, 0, 0/)
          planarityDistEq(1,1:3) = (/ 5.2, 3.7, 2.7/)
!         A
          planarityDistEq(2,1:3) = (/ 4.5, 4.0, 2.8/)
!         C
          planarityDistEq(3,1:3) = (/ 1.9, 0.79, 0.38/)
!         U
          planarityDistEq(4,1:3) = (/ 2.7, 1.4, 0.30/)
          
      alpam = alpam+pi
      alpbm = alpbm+pi
      calpam = cos(alpam)
      salpam = sin(alpam)
      calpbm = cos(alpbm)
      salpbm = sin(alpbm)

   END SUBROUTINE FILL_DNA_HB_PARAMS  
    
END MODULE DNA_HB_PARAMS
