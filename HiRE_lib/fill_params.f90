!> @file
!> Contains FILL_PARAMS module handling nucleic acid parameters
!>
!> Module to deal with NA paramters
!> @brief
!>
!> Contains the function to assign all nucleic acid parameters, as well as the required allocation and deallocation routines.
MODULE FILL_PARAMS
   USE PREC_HIRE
   IMPLICIT NONE

   CONTAINS
      !> Subroutine to fill nucleic acid parameters
      !> @brief
      !>
      !> This routine assigns the NA relevant from the scale.dat input.\n
      !> We then allocate the associated arrays and assign the nucleobase type informations 
      !> based on the residue information from the topology.\n
      !> Then the parameters are populated for RNA and DNA.
      !>
      !> @see ALLOC_NAPARAMS
      !> @see FILL_RNA_HB_PARAMS
      !> @see FILL_RNA_HB_PARAMS
      SUBROUTINE FILL_HIRE_PARAMS()
         USE VAR_DEFS, ONLY: NRES, RESFINAL, RESNAMES, RESTYPES
         USE UTILS_IO, ONLY: GETUNIT
         USE NAparams, ONLY: WC, WCCanonic, noWC, TIT, noWCq, SCORE_RNA, &
                           BLIST, BTYPE, BPROT, BP_CURR, &
                           cwwAA, cwwAG, cwwAC, cwwAU, cwwGC, cwwGU, cwwCC, cwwCU, cwwUU, &
                           !twwAA, twwAC, twwAU, twwGG, twwGU, twwGC, twwCC, twwCU, twwUU, &
                           cwh, twh, cws, tws, chh, thh, chs, ths, css, tss, tww
         USE RNA_HB_PARAMS, ONLY: FILL_RNA_HB_PARAMS
         USE DNA_HB_PARAMS, ONLY: FILL_DNA_HB_PARAMS
       
         IMPLICIT NONE
         INTEGER :: I !, BINFOUNIT
       
         !set parameters from SCORE
         !WC = SCORE_RNA(39)
         !WCCanonic = SCORE_RNA(40)
         !noWC = SCORE_RNA(41)



         TIT = SCORE_RNA(78)
         noWCq = SCORE_RNA(79)

         !twwAA = 2.00 * SCORE_RNA(61)
         !twwAC = 2.00 * SCORE_RNA(61)
         !twwAU = 2.10 * SCORE_RNA(61)
         !twwGG = 2.00 * SCORE_RNA(61)
         !twwGC = 1.80 * SCORE_RNA(61)
         !twwGU = 2.00 * SCORE_RNA(61)
         !twwCC = 2.40 * SCORE_RNA(61)
         !twwCU = 2.20 * SCORE_RNA(61)
         !twwUU = 2.30 * SCORE_RNA(61)

         !twwAA = 2.00 * SCORE_RNA(61)
         !twwAC = 2.00 * SCORE_RNA(61)
         !twwAU = 2.00 * SCORE_RNA(61)
         !twwGG = 2.00 * SCORE_RNA(61)
         !twwGC = 2.00 * SCORE_RNA(61)
         !twwGU = 2.00 * SCORE_RNA(61)
         !twwCC = 2.00 * SCORE_RNA(61)
         !twwCU = 2.00 * SCORE_RNA(61)
         !twwUU = 2.00 * SCORE_RNA(61)
         
         tww = SCORE_RNA(80)
         
         cwwAA = SCORE_RNA(81)
         cwwAG = SCORE_RNA(82)
         cwwAC = SCORE_RNA(83)
         cwwAU = SCORE_RNA(84)
         cwwGC = SCORE_RNA(85)
         cwwGU = SCORE_RNA(86)
         cwwCC = SCORE_RNA(87)
         cwwCU = SCORE_RNA(88)
         cwwUU = SCORE_RNA(89)
                  
         cwh = SCORE_RNA(90)
         twh = SCORE_RNA(91)
         cws = SCORE_RNA(92)
         tws = SCORE_RNA(93)
         chh = SCORE_RNA(94)
         thh = SCORE_RNA(95)
         chs = SCORE_RNA(96)
         ths = SCORE_RNA(97)
         css = SCORE_RNA(98)
         tss = SCORE_RNA(99)

         !fill base information
         CALL ALLOC_NAPARAMS(NRES)
         !QUERY: Should be able to add this to/derive it from the topology.
         !I think this information is now in the topology, but needs to be read out correctly.
         !BINFOUNIT = GETUNIT()
         !OPEN(UNIT=BINFOUNIT,FILE="bblist.dat",STATUS="old")   
         !READ(BINFOUNIT,*) (BLIST(I), BTYPE(I), BPROT(I), I =1, NRES)      
         !CLOSE(BINFOUNIT) 
         !The assignment of base type is hard coded - if these assignment are changed the base pair information is incorrect
         !DO NOT CHANGE UNLESS YOU CHANGE THE NA PARAMETER ASSIGNMENT!!!!!!!!
         BLIST(1:NRES) = RESFINAL(1:NRES)
         DO I=1,NRES
            IF ((RESNAMES(I)(1:1).EQ."G").OR.(RESNAMES(I)(1:2).EQ."DG")) THEN
               BTYPE(I) = 1
            ELSEIF ((RESNAMES(I)(1:1).EQ."C").OR.(RESNAMES(I)(1:2).EQ."DC")) THEN 
               BTYPE(I) = 3
            ELSEIF ((RESNAMES(I)(1:1).EQ."A").OR.(RESNAMES(I)(1:2).EQ."DA")) THEN 
               BTYPE(I) = 2
            ELSEIF ((RESNAMES(I)(1:1).EQ."U").OR.(RESNAMES(I)(1:2).EQ."DT")) THEN
               BTYPE(I) = 4
            ENDIF
            ! set residue types (RNA or DNA)
            IF (RESNAMES(I)(1:1).EQ."D") THEN
               RESTYPES(I) = 1
            ELSE 
               RESTYPES(I) = 0
            ENDIF
         ENDDO
         !no protonation for now
         BPROT(1:NRES) = 0

         !call parameter filling for RNA
         CALL FILL_RNA_HB_PARAMS()
         !call parameter filling for DNA
         CALL FILL_DNA_HB_PARAMS()
         BP_CURR(:,:) = .FALSE.
      END SUBROUTINE FILL_HIRE_PARAMS
    
      !> Allocate NA parameter arrays
      !>
      !> @param[in] NSIZE - number of nucleotides
      SUBROUTINE ALLOC_NAPARAMS(NSIZE)
         USE NAparams, ONLY: BLIST, BTYPE, BPROT, BOCC, BPCH, BP_CURR
         IMPLICIT NONE
         INTEGER, INTENT(IN) ::  NSIZE
         
         CALL DEALLOC_NAPARAMS()
         ALLOCATE(BLIST(NSIZE),BTYPE(NSIZE),BPROT(NSIZE),BOCC(NSIZE), &
                  BPCH(NSIZE), BP_CURR(NSIZE,NSIZE))
      END SUBROUTINE ALLOC_NAPARAMS 

      !> Deallocate NA parameter arrays 
      !> @brief   
      SUBROUTINE DEALLOC_NAPARAMS()
         USE NAparams, ONLY: BLIST, BTYPE, BPROT, BOCC, BPCH, BP_CURR    
         IF (ALLOCATED(BLIST)) DEALLOCATE(BLIST)
         IF (ALLOCATED(BTYPE)) DEALLOCATE(BTYPE)
         IF (ALLOCATED(BPROT)) DEALLOCATE(BPROT)
         IF (ALLOCATED(BOCC)) DEALLOCATE(BOCC)
         IF (ALLOCATED(BPCH)) DEALLOCATE(BPCH)
         IF (ALLOCATED(BP_CURR)) DEALLOCATE(BP_CURR) 
      END SUBROUTINE DEALLOC_NAPARAMS
      
END MODULE FILL_PARAMS


