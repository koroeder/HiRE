!> @file
!> Contains the modules NAparams, RNA_HB_PARAMS and DNA_HB_PARAMS

!> Module containing the global variables needed to describe RNA and DNA
MODULE NAparams
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE

   !> number of parameters given in scale.dat
   INTEGER, PARAMETER :: SCORESIZE=99
   !> HiRE potential parameters given in scale.dat
   REAL(KIND = REAL64) :: SCORE_RNA(SCORESIZE)
   
   !RNA bases
   !> Index of the final particle for each base
   INTEGER, ALLOCATABLE :: BLIST(:)
   !> Type of nulceobase
   INTEGER, ALLOCATABLE :: BTYPE(:)
   !> Protonation state of bases
   INTEGER, ALLOCATABLE :: BPROT(:)   
   !> Occupation state of bases
   INTEGER, ALLOCATABLE :: BOCC(:)
   !> QUERY - what is this?
   REAL(KIND = REAL64), ALLOCATABLE :: BPCH(:)
   !> Matrix of base pairing
   LOGICAL, ALLOCATABLE :: BP_CURR(:,:)
   
   !Cut off definitions
   !> NO CUTOFF
   REAL(KIND = REAL64), PARAMETER :: CUT   = 1.0D2
   !> divide the 1-4 VDW interactions by 8.0
   REAL(KIND = REAL64), PARAMETER :: SCNB  = 8.0D1 
   !> divide the 1-4 ELEC interactions by 2.0 - QUERY this is the same as SCNB at the moment?!?
   REAL(KIND = REAL64), PARAMETER :: SCEE  = 8.0D1 
   !> 2.0 dielectric constant epsilon = 2r - QUERY it is zero?!?
   REAL(KIND = REAL64), PARAMETER :: IDIEL = 0.0D0 
   !> Dielectric constant
   REAL(KIND = REAL64), PARAMETER :: DIELC = 1.0D0
   !> Zero to be used in filling parameter arrays         
   REAL(KIND = REAL64), PARAMETER :: Z = 0.0D0
   
   !> Possible base names
   CHARACTER(LEN=1), DIMENSION(4) :: basenum = (/'G','A','C','U'/)        !
   
   !> Parameter array
   INTEGER :: NPARAM(4,4) !
   
   !> rna/dna parameter: wc coef
   REAL(KIND = REAL64) :: wc                 
   !> rna/dna parameter: wcCanonic coef
   REAL(KIND = REAL64) :: wcCanonic           
   !> rna/dna parameter: noWc coef
   REAL(KIND = REAL64) :: noWc               
   !> rna parameter titration
   REAL(KIND = REAL64) :: tit                
   !> rna parameter hb charged base
   REAL(KIND = REAL64) :: noWCq 
   !> cis Watson-Crick base pairs for  RNA            
   REAL(KIND = REAL64) :: cwwAA
   REAL(KIND = REAL64) :: cwwAG 
   REAL(KIND = REAL64) :: cwwAC 
   REAL(KIND = REAL64) :: cwwAU 
   REAL(KIND = REAL64) :: cwwGC 
   REAL(KIND = REAL64) :: cwwGU 
   REAL(KIND = REAL64) :: cwwCC 
   REAL(KIND = REAL64) :: cwwCU 
   REAL(KIND = REAL64) :: cwwUU
   !> non-canonical base pairs for RNA:
   !> c/t - cis/trans
   !> w - Watson-Crick, h - Hoogsteen, s - Sugar edge
   REAL(KIND = REAL64) :: cwh 
   REAL(KIND = REAL64) :: twh 
   REAL(KIND = REAL64) :: cws 
   REAL(KIND = REAL64) :: tws 
   REAL(KIND = REAL64) :: chh 
   REAL(KIND = REAL64) :: thh 
   REAL(KIND = REAL64) :: chs 
   REAL(KIND = REAL64) :: ths
   REAL(KIND = REAL64) :: css 
   REAL(KIND = REAL64) :: tss
   REAL(KIND = REAL64) :: tww
      
   !Variables for titration
   !> titration start
   INTEGER, ALLOCATABLE :: fatortit(:)
   
   !RNA cutoffs - QUERY - what are all these cutoffs?
   !> What is this?  
   REAL(KIND = REAL64), PARAMETER :: rcut2_caca_scsc_out = 1.44D2
   !> What is this?
   REAL(KIND = REAL64), PARAMETER :: rcut2_caca_scsc_in = 1.00D2
   !> Hydrogen bonding interactions cutoff
   REAL(KIND = REAL64), PARAMETER :: rcut2_hb_mcmc_out = 1.44D2
   !> What is this?
   REAL(KIND = REAL64), PARAMETER :: rcut2_hb_mcmc_in = 1.00D2
   !> What is this?
   REAL(KIND = REAL64), PARAMETER :: rcut2_4b_out = 1.44D2
   !> What is this?
   REAL(KIND = REAL64), PARAMETER :: rcut2_4b_in = 1.00D2
   !> What is this?
   REAL(KIND = REAL64), PARAMETER :: rcut2_lj_out = 1.44D2
   !> What is this?
   REAL(KIND = REAL64), PARAMETER :: rcut2_lj_in = 1.00D2
   
   !Additional cutoffs used for speed of computation
   !> Cut off for DH calculations
   REAL(KIND = REAL64), PARAMETER :: DHCUT = 35.0
   !> Cut off for non-bonded interactions (square of distance)
   REAL(KIND = REAL64), PARAMETER :: NBCUT = 400.0
   SAVE
END MODULE NAparams

!> Module to fill the RNA paramters
MODULE RNA_HB_PARAMS
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE
   !> Parameter array
   INTEGER :: NPARAM(4,4) !
   !> not used
   REAL(KIND = REAL64) :: alpa(4,4)            ! 
   !> not used
   REAL(KIND = REAL64) :: alpb(4,4)            !
   !> Parameters for base A
   REAL(KIND = REAL64) :: alpam(15,4,4)         ! 
   !> Parameters for base B
   REAL(KIND = REAL64) :: alpbm(15,4,4)         !
   !> Scaling for cosine term for angle in base A
   REAL(KIND = REAL64) :: calpam(15,4,4)        ! 
   !> Scaling for cosine term for angle in base B
   REAL(KIND = REAL64) :: calpbm(15,4,4)        !
   !> Scaling for sine term for angle in base A
   REAL(KIND = REAL64) :: salpam(15,4,4)        ! 
   !> Scaling for sine term for angle in base B
   REAL(KIND = REAL64) :: salpbm(15,4,4)        !
   !> Reference value for h-bonding distance
   REAL(KIND = REAL64) :: dREF(15,4,4)          ! 

   !> All WC and noWC parameters ready for look up
   REAL(KIND = REAL64) :: s(15,4,4,2,2)         !
   !> Equilibrium distances for planarity term
   REAL(KIND = REAL64) :: planarityDistEq(4,3) !

   CONTAINS
      
      !> Routine to fill module arrays for RNA
      SUBROUTINE FILL_RNA_HB_PARAMS
         USE NAparams, ONLY: WC, WCCanonic, noWC, TIT, noWCq, Z, &
            cwwAA, cwwAG, cwwAC, cwwAU, cwwGC, cwwGU, cwwCC, cwwCU, cwwUU, &
            cwh, twh, cws, tws, chh, thh, chs, ths, css, tss, tww
       
         USE NUM_DEFS, ONLY: PI
         
         ! non WW bp scales are multipied by the number of HB formed 
         
         alpam(:,:,:) = 2.0D0
         alpbm(:,:,:) = 2.0D0
         !         A-A
         dREF(1:9,2,2) =  (/ 6.585, 6.618, 6.956, 5.424, 6.664, 6.661, 6.006, 5.917, 5.389/)       
			   !cHS_AA, cSS_AA, cWS_AA, cWW_Aa, tHS_AA, tSS_AA, tWH_AA, tWS_AA, tWW_AA                    
         alpam(1:9,2,2) = (/ -1.654, 1.494, -1.095, 0.112, -1.912, 0.768, -0.300, -0.473, -0.568 /)
         alpbm(1:9,2,2) = (/ 1.531, 0.776, 1.595, -0.647, 1.400, 1.679, -1.889, 1.212, -0.566 /)
         s(1:9,2,2,1,1) = (/0.5*chs, 1.5*css, cws, cwwAA, 1.5*ths, 1.5*tss, twh, tws, tww/)    	!Br2  1,1 -> q1=0, q2=0
         s(1:9,2,2,1,2) = (/0.5*chs, 1.5*css, cws, z,     1.5*ths, 1.5*tss, twh, tws, z/)       !Br2  1,2 -> q1=0, q2=1
         s(1:9,2,2,2,1) = (/0.5*chs, 1.5*css, z,   z,     1.5*ths, z,       z,	 z ,  z/)  	!Br2  2,1 -> q1=1, q2=0
         s(1:9,2,2,2,2) = (/0.5*chs, 1.5*css, z,   z,     1.5*ths, z,       z,	 z,   z /)        !Br2  2,2 -> q1=1, q2=1
         Nparam(2,2) = 9  !Br3  CHECK angle sign when order is inverted when one is negative

         !         A-C
         dREF(1:14,2,3) = (/6.341, 8.801, 6.526, 7.184, 6.794, 6.740, 5.589, 6.883, 7.039, 6.776, 6.845, 5.813, 7.064, 5.281 /)         
			! cHS_AC, cSS_AC, tHH_AC, tHS_AC, tSS_AC, tWS_AC, tWW_AC, cHS_CA, cSS_CA, cWS_CA, tHS_CA, tWH_CA, tWS_CA, cWW_CA
         alpam(1:14,2,3) = (/-1.858, 1.727, -2.062, -1.719, 0.695, -0.758, -0.728, 1.680, 0.998, 1.631, 1.372, -1.794, 1.536, -0.576/)
         alpbm(1:14,2,3) = (/-2.363, -1.758, 1.228, -2.186, -2.462, -2.142, -0.387, 0.566, -2.272, -0.262, 0.660, -0.815, -0.295, -1.321/)   
         s(1:14,2,3,1,1) = (/ 0.5*chs, css, 0.5*thh, ths, tss, tws, tww, 0.5*chs, css, cws, ths, twh, tws, cwwAC /)    	!Br2  1,1 -> q1=0, q2=0
         s(1:14,2,3,1,2) = (/ 0.5*chs, css, 0.5*thh, ths, tss, tws, z,   0.5*chs, css, cws, ths, twh, tws, z   /)         !Br2  1,2 -> q1=0, q2=1
         s(1:14,2,3,2,1) = (/ 0.5*chs, css, 0.5*thh, ths, z,   tws, z,   0.5*chs, css, z,   ths, z,   tws, z   /)      	!Br2  2,1 -> q1=1, q2=0
         s(1:14,2,3,2,2) = (/ 0.5*chs, css, 0.5*thh, ths, z,   tws, z,   0.5*chs, css, z,   ths, z,   tws, z   /)         !Br2  2,2 -> q1=1, q2=1

         Nparam(2,3) = 14
         dREF(:,3,2) = dREF(:,2,3)
         alpam(:,3,2) = alpbm(:,2,3)
         alpbm(:,3,2) = alpam(:,2,3)
         s(:,3,2,1,1) = s(:,2,3,1,1)
         s(:,3,2,1,2) = s(:,2,3,2,1)
         s(:,3,2,2,1) = s(:,2,3,1,2)
         s(:,3,2,2,2) = s(:,2,3,2,2)
         Nparam(3,2) = Nparam(2,3)

         !         A-G
         dREF(1:15,2,1) = (/ 6.980, 7.042, 7.035, 6.049, 7.588, 6.840, 6.610, 6.147, 5.979, 7.340, 6.560, 6.042, 5.375, 6.519, 4.773 /)         
			!  cHH_AG,  cHS_AG, cSS_AG, cWH_AG, cWS_AG, tHH_AG, tHS_AG, tSS_AG, tWH_AG, cHS_GA, cSS_GA, cWH_GA, cWS_GA, tSS_GA, cWW_GA
         alpam(1:15,2,1) = (/ -1.644, -2.231, 1.539, -0.506, -1.166, -1.969, -1.873, 1.099, -0.173, 1.092, 0.662, -2.226, 1.304, 1.361, -0.092/)
         alpbm(1:15,2,1) = (/ -1.190, 0.810, 0.837, -2.247, 1.667, -1.714, 1.476, 1.537, -1.968, -1.754, 1.501, -0.860, -0.432, 1.410, -0.444/)
         s(1:15,2,1,1,1) = (/ 0.5*chh, 0.5*chs, 1.5*css, cwh, cws, 0.5*thh, ths, 1.5*tss, twh, 0.5*chs, 2*css, cwh, 1.5*cws, 1.5*tss, cwwAG /)       	!Br2  1,1 -> q1=0, q2=0
         s(1:15,2,1,1,2) = (/ 0.5*chh, 0.5*chs, 1.5*css, cwh, cws, 0.5*thh, ths, 1.5*tss, twh, 0.5*chs, z,     cwh, 1.5*cws, 1.5*tss, z /)        !Br2  1,2 -> q1=0, q2=1
         s(1:15,2,1,2,1) = (/ 0.5*chh, 0.5*chs, 1.5*css, z,   z,   0.5*thh, ths, z,       z,   0.5*chs, 2*css, z,   z,       1.5*tss, z /)    		!Br2  2,1 -> q1=1, q2=0
         s(1:15,2,1,2,2) = (/ 0.5*chh, 0.5*chs, 1.5*css, z,   z,   0.5*thh, ths, z,       z,   0.5*chs, z,     z,   z,       1.5*tss, z /)    		!Br2  2,2 -> q1=1, q2=1

         Nparam(2,1) = 15
         dREF(:,1,2) = dREF(:,2,1)
         alpam(:,1,2) = alpbm(:,2,1)
         alpbm(:,1,2) = alpam(:,2,1)
         s(:,1,2,1,1) = s(:,2,1,1,1)
         s(:,1,2,1,2) = s(:,2,1,2,1)
         s(:,1,2,2,1) = s(:,2,1,1,2)
         s(:,1,2,2,2) = s(:,2,1,2,2)
         Nparam(1,2) = Nparam(2,1)

         !         A-U
         dREF(1:15,2,4) = (/ 7.416, 8.684, 6.939, 4.807, 7.253, 7.186, 6.812, 4.765, 6.824, 7.037, 5.941, 5.670, 7.161, 5.888, 5.613 /)        
			  ! cHS_AU, cSS_AU, cWS_AU, cWW_AU, tHH_AU, tHS_AU, tSS_AU, tWW_AU, cHS_UA, cSS_UA, cWH_UA, cWS_UA, tHS_UA, tWH_UA, tWS_UA
         alpam(1:15,2,4) = (/ -1.591, 1.584, -0.726, -0.140, -1.734, -1.687, 0.765, -0.190, 1.308, 1.042, -2.229, 1.503, 0.848, -2.196, 1.512/)
         alpbm(1:15,2,4) = (/ -2.165, -1.779, -2.397, -0.911,  0.471, -2.215, -2.448, -0.937, 0.100, -2.281, -0.697, -0.619, 0.302, -1.264,-1.048 /)
         s(1:15,2,4,1,1) = (/0.5*chs, css, cws, cwwAU, 0.5*thh, ths, tss, 1.5*tww, 0.5*chs, 1.5*css, 1.5*cwh, 1.5*cws, 0.5*ths, twh, tws/) 	!Br2  1,1 -> q1=0, q2=0 
                                             												!OKKIO : cambiato a mano 14.4 -> 16!!!
         s(1:15,2,4,1,2) = (/0.5*chs, css, cws, z,     0.5*thh, ths, tss, z,       0.5*chs, 1.5*css, 1.5*cwh, 1.5*cws, 0.5*ths, twh, tws/)         !Br2  1,2 -> q1=0, q2=1
         s(1:15,2,4,2,1) = (/0.5*chs, css, z,   z,     0.5*thh, ths, z,   z,       0.5*chs, 1.5*css, z,       z,       0.5*ths, z,   z/)           !Br2  2,1 -> q1=1, q2=0
         s(1:15,2,4,2,2) = (/0.5*chs, css, z,   z,     0.5*thh, ths, z,   z,       0.5*chs, 1.5*css, z,       z,       0.5*ths, z,   z/)           !Br2  2,2 -> q1=1, q2=1
         Nparam(2,4) = 15
            
         dREF(:,4,2) = dREF(:,2,4)
         alpam(:,4,2) = alpbm(:,2,4)
         alpbm(:,4,2) = alpam(:,2,4)
         s(:,4,2,1,1) = s(:,2,4,1,1)
         s(:,4,2,1,2) = s(:,2,4,2,1)
         s(:,4,2,2,1) = s(:,2,4,1,2)
         s(:,4,2,2,2) = s(:,2,4,2,2)
         Nparam(4,2) = Nparam(2,4)

        !       C-C
        dREF(1:9,3,3) =  (/ 6.892, 8.874, 6.136, 7.046, 5.350, 7.224, 6.162, 7.260, 4.744 /)  
        !   cHS_CC, cSS_CC, cWH_CC, cWS_CC, cWW_cC, tHS_CC, tWH_CC, tWS_CC, tWW_CC
        alpam(1:9,3,3) = (/ 1.015, -1.855, -1.015, -0.417, -0.442, 0.639, -1.114, -0.336, 0.903 /)
        alpbm(1:9,3,3) = (/ -2.226, -2.592, 0.738, -2.459, -1.554, -2.150, 1.026, -1.886, -0.939  /)
        s(1:9,3,3,1,1) = (/ chs, css, cwh, cws, cwwCC, ths, twh, tws, tww /) 		!Br2  1,1 -> q1=0, q2=0
        s(1:9,3,3,1,2) = (/ chs, css, cwh, cws, z,     ths, twh, tws, tww  /) 		!Br2  1,2 -> q1=0, q2=1
        s(1:9,3,3,2,1) = (/ chs, css, z,   z,   z,     ths, z,   tws, tww  /) 		!Br2  2,1 -> q1=1, q2=0
        s(1:9,3,3,2,2) = (/ chs, css, z,   z,   z,     ths, z,   tws, tww /) 		!Br2  2,2 -> q1=1, q2=1
        Nparam(3,3) = 9

         !         C-G
         dREF(1:14,3,1) = (/  8.206, 7.563, 6.419, 6.779, 5.981, 5.762, 7.218, 10.351, 5.800, 4.801, 7.177, 7.517, 5.710, 5.231 /)  
			 !  cHS_CG, cSS_CG, cWH_CG, cWS_CG, tWH_CG, tWS_CG, cHH_GC, cSS_GC, cWS_GC, cWW_GC, tHH_GC, tSS_GC, tWS_GC , tWW_GC 
         alpam(1:14,3,1) = (/  0.131, -2.343, -0.620, -0.231, -0.786, -0.450, -2.354, -2.053, -1.911, -0.944, 0.991, -2.076, -2.059, -1.489 /)
         alpbm(1:14,3,1) = (/  1.645, 1.014, -1.914, 1.639, -1.861, 1.207, -2.054, 1.759, -0.408, -0.388, -2.207, 1.074, -0.220, -0.012 /)
         s(1:14,3,1,1,1) = (/ 0.5*chs, 1.5*css, cwh, cws, twh, 1.5*tws, 0.5*chh, css, cws, 1.5*cwwGC, 0.5*thh, tss, tws, tww/)       		!Br2  1,1 -> q1=0, q2=0
         s(1:14,3,1,1,2) = (/ 0.5*chs, 1.5*css, cwh, cws, twh, 1.5*tws, 0.5*chh, css, cws, z,     0.5*thh, tss, tws, z /)                  !Br2  1,2 -> q1=0, q2=1
         s(1:14,3,1,2,1) = (/ 0.5*chs, 1.5*css, z,   z,   z,   z,       0.5*chh, css, z,   z,     0.5*thh, tss, z,   z  /)                  	!Br2  2,1 -> q1=1, q2=0
         s(1:14,3,1,2,2) = (/ 0.5*chs, 1.5*css, z,   z,   z,   z,       0.5*chh, css, z,   z,     0.5*thh, tss, z,   z  /)                   !Br2  2,2 -> q1=1, q2=1
         Nparam(3,1) = 14
            
         dREF(:,1,3) = dREF(:,3,1)
         alpam(:,1,3) = alpbm(:,3,1)
         alpbm(:,1,3) = alpam(:,3,1)
         s(:,1,3,1,1) = s(:,3,1,1,1)
         s(:,1,3,1,2) = s(:,3,1,2,1)
         s(:,1,3,2,1) = s(:,3,1,1,2)
         s(:,1,3,2,2) = s(:,3,1,2,2)
         Nparam(1,3) = Nparam(3,1)

         !         C-U
         dREF(1:9,3,4) = (/ 7.232, 9.669, 7.027, 6.363, 6.859, 7.406, 4.899, 7.568, 5.931 /)             
			 ! cHS_CU, cSS_CU, cWS_CU, cWW_CU, tHH_CU, tHS_CU, tWW_CU, cSS_UC, tWS_UC
         alpam(1:9,3,4) = (/0.779, -2.662, -0.392, -0.550, 0.559, 0.635, -0.901, -1.624, -1.904 /)
         alpbm(1:9,3,4) = (/-2.587, -2.127, -2.430, -0.393, 0.481, -2.035, -0.924, -2.363, -0.793 /)            !Br3 CHECK when particles are inverted
         s(1:9,3,4,1,1) = (/0.5*chs, css, cws, cwwCU, 0.5*thh, ths, tww, css, tws/)              		!Br2  1,1 -> q1=0, q2=0     
         s(1:9,3,4,1,2) = (/0.5*chs, css, cws, cwwCU, 0.5*thh, ths, z,   css, tws/)               		!Br2  1,2 -> q1=0, q2=1
         s(1:9,3,4,2,1) = (/0.5*chs, css, z,   cwwCU, 0.5*thh, ths, z,   css, z /)               		!Br2  2,1 -> q1=1, q2=0
         s(1:9,3,4,2,2) = (/0.5*chs, css, z,   cwwCU, 0.5*thh, ths, z,   css, z /)               		!Br2  2,2 -> q1=1, q2=1
         Nparam(3,4) = 9
            
         dREF(:,4,3) = dREF(:,3,4)
         alpam(:,4,3) = alpbm(:,3,4)
         alpbm(:,4,3) = alpam(:,3,4)
         s(:,4,3,1,1) = s(:,3,4,1,1)
         s(:,4,3,1,2) = s(:,3,4,2,1)
         s(:,4,3,2,1) = s(:,3,4,1,2)
         s(:,4,3,2,2) = s(:,3,4,2,2)
         Nparam(4,3) = Nparam(3,4)

         !         G-G
         dREF(1:9,1,1) =  (/ 8.796, 7.240, 6.990, 6.661, 9.127, 7.333, 6.135, 6.251, 5.176 /)        
			  ! cHH_gG, cHS_GG, cSS_GG, cWS_GG, tHH_gG, tHS_GG, tSS_GG, tWH_GG, tWW_GG
         alpam(1:9,1,1) = (/ 2.979, -1.883, 1.575, -0.851, -3.031, -1.821, 1.330, -0.120, -0.881 /)
         alpbm(1:9,1,1) = (/ -1.953, 1.062, 0.920, 1.556, -1.773, 1.039, 1.328, -2.322, -0.871 /)
         s(1:9,1,1,1,1) = (/ 0.5*chh, 0.5*chs, 1.5*css, cws, 0.5*thh, 0.5*ths, 2*tss, twh, tww /)      		!Br2  1,1 -> q1=0, q2=0     
         s(1:9,1,1,1,2) = (/ 0.5*chh, 0.5*chs, 1.5*css, cws, 0.5*thh, 0.5*ths, 2*tss, twh, z /)               	!Br2  1,2 -> q1=0, q2=1
         s(1:9,1,1,2,1) = (/ 0.5*chh, 0.5*chs, 1.5*css, cws, 0.5*thh, 0.5*ths, 2*tss, z,   z /)   			!Br2  2,1 -> q1=1, q2=0
         s(1:9,1,1,2,2) = (/ 0.5*chh, 0.5*chs, 1.5*css, cws, 0.5*thh, 0.5*ths, 2*tss, z,   z /)                	!Br2  2,2 -> q1=1, q2=1
         Nparam(1,1) = 9

         !         G-U
         dREF(1:14,1,4) = (/ 8.955, 5.826, 7.614, 5.739, 5.756, 5.413, 6.677, 7.537, 6.764, 5.544, 6.985, 7.162, 6.873,  5.263 /)              
			! cSS_GU, cWS_GU, tSS_GU, tWH_GU, tWS_GU, tWW_GU, cHS_UG, cSS_UG, cWH_UG, cWS_UG, tHS_UG, tWH_UG, tWS_UG, cWW_UG   
         alpam(1:14,1,4) = (/ 1.761, -0.380, 0.949, -0.465, -0.333, -0.964, 1.476, 1.025, -2.531, 1.590, 1.151, -2.475, 1.293, -0.779  /)
         alpbm(1:14,1,4) = (/ -1.630, -1.824, -2.298, 0.614, -2.159, -0.442, 0.676, -2.322, -0.829, -0.650, 0.452, -0.875, -1.623, -1.470  /) 
         s(1:14,1,4,1,1) = (/ css, cws, tss, twh, tws, tww, 0.5*chs, 1.5*css, cwh, cws, 0.5*ths, twh, 0.5*tws, cwwGU /)     		!Br2  1,1 -> q1=0, q2=0 !OKKIO : cambiato a mano 14.7 -> 16!!!
         s(1:14,1,4,1,2) = (/ css, cws, tss, z,   tws, tww, 0.5*chs, 1.5*css, cwh, cws, 0.5*ths, twh, 0.5*tws, z /)                 	!Br2  1,2 -> q1=0, q2=1
         s(1:14,1,4,2,1) = (/ css, z,   tss, z,   z,   tww, 0.5*chs, 1.5*css, z,   z,   0.5*ths, z,   0.5*tws, z  /)                 	!Br2  2,1 -> q1=1, q2=0
         s(1:14,1,4,2,2) = (/ css, z,   tss, z,   z,   tww, 0.5*chs, 1.5*css, z,   z,   0.5*ths, z,   0.5*tws, z /)                 	!Br2  2,2 -> q1=1, q2=1
         Nparam(1,4) = 14
            
         dREF(:,4,1) = dREF(:,1,4)
         alpam(:,4,1) = alpbm(:,1,4)
         alpbm(:,4,1) = alpam(:,1,4)
         s(:,4,1,1,1) = s(:,1,4,1,1)
         s(:,4,1,1,2) = s(:,1,4,2,1)
         s(:,4,1,2,1) = s(:,1,4,1,2)
         s(:,4,1,2,2) = s(:,1,4,2,2)
         Nparam(4,1) = Nparam(1,4)

         !         U-U
         dREF(1:7,4,4) =  (/ 7.177, 8.091, 6.311, 5.306, 5.898, 6.364 , 5.174 /)           
			  ! cHS_UU, cSS_UU, cWS_UU, cWW_uU, tWH_UU, tWS_UU, tWW_UU
         alpam(1:7,4,4) = (/ 1.151, -2.372, -1.079, -0.501, -0.806, -1.339, -0.438 /)
         alpbm(1:7,4,4) = (/ -2.380, -1.833, -2.169, -1.337, 0.473, -2.022, -0.516 /)
         s(1:7,4,4,1,1) = (/ 0.5*chs, css, cws, cwwUU, twh, tws, tww/)        		!Br2  1,1 -> q1=0, q2=0
         s(1:7,4,4,1,2) = (/ 0.5*chs, css, cws, z,     twh, tws, z/)                	!Br2  1,2 -> q1=0, q2=1
         s(1:7,4,4,2,1) = (/ 0.5*chs, css, z,   z,     z,   z,    z/)                    	!Br2  2,1 -> q1=1, q2=0
         s(1:7,4,4,2,2) = (/ 0.5*chs, css, z,   z,     z,   z,    z/)                    	!Br2  2,2 -> q1=1, q2=1
         Nparam(4,4) = 7


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
!         alpam = alpam+pi
!         alpbm = alpbm+pi
         calpam = cos(alpam)
         salpam = sin(alpam)
         calpbm = cos(alpbm)
         salpbm = sin(alpbm)  
      END SUBROUTINE FILL_RNA_HB_PARAMS
END MODULE RNA_HB_PARAMS

!> Module to provide DNA parameters
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
            s(1:3,2,4,1,1) = (/ 2.2*wcCanonic, 2*noWc, 2*noWc/)   !Br2  1,1 -> q1=0, q2=0   !OKKIO : cambiato a mano 14.4 -> 16!
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
