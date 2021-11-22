MODULE FILL_PARAMS
  USE PREC_HIRE
  IMPLICIT NONE

  CONTAINS
  
    SUBROUTINE FILL_HIRE_PARAMS()
       USE VAR_DEFS, ONLY: NRES
       USE UTILS_IO, ONLY: GETUNIT
       USE NAparams, ONLY: WC, WCCanonic, noWC, TIT, noWCq, SCORE_RNA, &
                           BLIST, BTYPE, BPROT, BP_CURR
       USE RNA_HB_PARAMS, ONLY: FILL_RNA_HB_PARAMS
       USE DNA_HB_PARAMS, ONLY: FILL_DNA_HB_PARAMS
       
       IMPLICIT NONE
       INTEGER :: BINFOUNIT, I
       
       !set parameters from SCORE
       WC = SCORE_RNA(39)
       WCCanonic = SCORE_RNA(40)
       noWC = SCORE_RNA(41)
       TIT = SCORE_RNA(42)
       noWCq = SCORE_RNA(43)
       
       !fill base information
       CALL ALLOC_NAPARAMS(NRES)
       !QUERY: Should be able to add this to/derive it from the topology.
       !I thionk this information is now in the topology, but needs to be read out correctly.
       BINFOUNIT = GETUNIT()
       OPEN(UNIT=BINFOUNIT,FILE="bblist.dat",STATUS="old")   
       READ(BINFOUNIT,*) (BLIST(I), BTYPE(I), BPROT(I), I =1, NRES)       
          
       !call parameter filling for RNA
       CALL FILL_RNA_HB_PARAMS()
       !call parameter filling for DNA
       CALL FILL_DNA_HB_PARAMS()
       BP_CURR(:,:) = .FALSE.
    END SUBROUTINE FILL_HIRE_PARAMS
    
    SUBROUTINE ALLOC_NAPARAMS(NSIZE)
       USE NAparams, ONLY: BLIST, BTYPE, BPROT, BOCC, BPCH, BP_CURR
       IMPLICIT NONE
       INTEGER, INTENT(IN) ::  NSIZE
       
       CALL DEALLOC_NAPARAMS()
       ALLOCATE(BLIST(NSIZE),BTYPE(NSIZE),BPROT(NSIZE),BOCC(NSIZE), &
                BPCH(NSIZE), BP_CURR(NSIZE,NSIZE))
    END SUBROUTINE ALLOC_NAPARAMS 

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


