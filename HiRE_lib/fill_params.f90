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
                           BLIST, BTYPE, BPROT, BP_CURR
         USE RNA_HB_PARAMS, ONLY: FILL_RNA_HB_PARAMS
         USE DNA_HB_PARAMS, ONLY: FILL_DNA_HB_PARAMS
       
         IMPLICIT NONE
         INTEGER :: I !, BINFOUNIT
       
         !set parameters from SCORE
         WC = SCORE_RNA(39)
         WCCanonic = SCORE_RNA(40)
         noWC = SCORE_RNA(41)
         TIT = SCORE_RNA(42)
         noWCq = SCORE_RNA(43)
         
         !fill base information
         CALL ALLOC_NAPARAMS(NRES)
         !QUERY: Should be able to add this to/derive it from the topology.
         !I think this information is now in the topology, but needs to be read out correctly.
         !BINFOUNIT = GETUNIT()
         !OPEN(UNIT=BINFOUNIT,FILE="bblist.dat",STATUS="old")   
         !READ(BINFOUNIT,*) (BLIST(I), BTYPE(I), BPROT(I), I =1, NRES)      
         !CLOSE(BINFOUNIT) 
         BLIST(1:NRES) = RESFINAL(1:NRES)
         DO I=1,NRES
            IF (RESNAMES(I)(1:1).EQ."G") THEN
               BTYPE(I) = 1
               RESTYPES(I) = 0
            ELSEIF (RESNAMES(I)(1:1).EQ."C") THEN 
               BTYPE(I) = 2
               RESTYPES(I) = 0
            ELSEIF (RESNAMES(I)(1:1).EQ."A") THEN 
               BTYPE(I) = 3
               RESTYPES(I) = 0
            ELSEIF (RESNAMES(I)(1:1).EQ."U") THEN
               BTYPE(I) = 4
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


