MODULE CG_DATA
   USE TOP_GLOBALS
   USE PREC_HIRE

   CONTAINS

      !get CG residue names and termini location
      SUBROUTINE GET_RES_DATA(SEQ)
         CHARACTER(LEN=4), INTENT(IN) :: SEQ(NRES)

         CHARACTER(LEN=4) :: CURRNAME
         INTEGER :: NDUMMYTER
         INTEGER :: DUMMYTER(NRES,2)
         INTEGER :: J
         INTEGER :: ATOMCOUNTER
         
         ALLOCATE(CGRESNAMES(NRES),RESTYPE(NRES),CGSTART(NRES),CGFINAL(NRES))
         NDUMMYTER = 0
         ATOMCOUNTER = 0
         RESTYPE(1:NRES) = 0 !default is RNA
         DO J=1,NRES
            CGSTART(J) = ATOMCOUNTER + 1
            CURRNAME = SEQ(J)
            !check whether residue name has a D in it (DNA)
            IF (INDEX(CURRNAME,"D").GT.0) THEN
               RESTYPE(J) = 1
            END IF
            !ion types
            IF (CURRNAME.EQ."MG") THEN
               CGRESNAMES(J) = "MG"
               ATOMCOUNTER = ATOMCOUNTER + 1
               RESTYPE(J) = 2
            ELSE IF (CURRNAME.EQ."NA") THEN
               CGRESNAMES(J) = "NA"
               ATOMCOUNTER = ATOMCOUNTER + 1
               RESTYPE(J) = 2
            ELSE IF (CURRNAME.EQ."K") THEN
               CGRESNAMES(J) = "K"
               ATOMCOUNTER = ATOMCOUNTER + 1
               RESTYPE(J) = 2
            ELSE IF (CURRNAME.EQ."CL") THEN
               CGRESNAMES(J) = "CL"
               ATOMCOUNTER = ATOMCOUNTER + 1  
               RESTYPE(J) = 2      
            !check which residue we have
            ELSE IF (INDEX(CURRNAME,"A").GT.0) THEN
               IF (INDEX(CURRNAME,"5").GT.0) THEN
                  CGRESNAMES(J) = "A5"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               ELSE IF (INDEX(CURRNAME,"3").GT.0) THEN                
                  CGRESNAMES(J) = "A3"
                  ATOMCOUNTER = ATOMCOUNTER + 7
               ELSE
                  CGRESNAMES(J) = "A"
                  ATOMCOUNTER = ATOMCOUNTER + 7
               END IF                 
            ELSE IF (INDEX(CURRNAME,"C").GT.0) THEN
               IF (INDEX(CURRNAME,"5").GT.0) THEN
                  CGRESNAMES(J) = "C5"
                  ATOMCOUNTER = ATOMCOUNTER + 5
               ELSE IF (INDEX(CURRNAME,"3").GT.0) THEN                
                  CGRESNAMES(J) = "C3"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               ELSE
                  CGRESNAMES(J) = "C"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               END IF  
            ELSE IF (INDEX(CURRNAME,"G").GT.0) THEN
               IF (INDEX(CURRNAME,"5").GT.0) THEN
                  CGRESNAMES(J) = "G5"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               ELSE IF (INDEX(CURRNAME,"3").GT.0) THEN                
                  CGRESNAMES(J) = "G3"
                  ATOMCOUNTER = ATOMCOUNTER + 7
               ELSE
                  CGRESNAMES(J) = "G"
                  ATOMCOUNTER = ATOMCOUNTER + 7
               END IF  
            ELSE IF (INDEX(CURRNAME,"T").GT.0) THEN
               IF (INDEX(CURRNAME,"5").GT.0) THEN
                  CGRESNAMES(J) = "T5"
                  ATOMCOUNTER = ATOMCOUNTER + 5
               ELSE IF (INDEX(CURRNAME,"3").GT.0) THEN                
                  CGRESNAMES(J) = "T3"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               ELSE
                  CGRESNAMES(J) = "T"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               END IF  
            ELSE IF (INDEX(CURRNAME,"U").GT.0) THEN
               IF (INDEX(CURRNAME,"5").GT.0) THEN
                  CGRESNAMES(J) = "U5"
                  ATOMCOUNTER = ATOMCOUNTER + 5
               ELSE IF (INDEX(CURRNAME,"3").GT.0) THEN                
                  CGRESNAMES(J) = "U3"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               ELSE
                  CGRESNAMES(J) = "U"
                  ATOMCOUNTER = ATOMCOUNTER + 6
               END IF  
            ELSE
               WRITE(*,*) "Residue name ", CURRNAME, " not recognised - STOP"
               STOP
            END IF
            !check for terminal position
            IF (INDEX(CURRNAME,"5").GT.0) THEN
               NDUMMYTER = NDUMMYTER + 1
               DUMMYTER(NDUMMYTER,1) = J
            ELSE IF (INDEX(CURRNAME,"3").GT.0) THEN
               DUMMYTER(NDUMMYTER,2) = J
            END IF
            CGFINAL(J) = ATOMCOUNTER
         END DO
         NTERMINI = NDUMMYTER
         ALLOCATE(TERMINI(NTERMINI,2))
         TERMINI(1:NTERMINI,1:2) = DUMMYTER(1:NTERMINI,1:2)
         NATOMS = ATOMCOUNTER
         ALLOCATE(CGNAMES(NATOMS),CGTYPE(NATOMS),CGMASS(NATOMS),CGCHARGE(NATOMS),XYZCG(3*NATOMS))
      END SUBROUTINE GET_RES_DATA

      SUBROUTINE ASSIGN_GRAIN_DATA()
         INTEGER :: J
         INTEGER :: ATOMCOUNTER

         ATOMCOUNTER = 0 
         DO J=1,NRES
            ! the type number is synced with the HiRE library
            IF (RESTYPE(J).EQ.2) THEN
               ATOMCOUNTER = ATOMCOUNTER + 1
               IF (CGRESNAMES(J).EQ."MG") THEN
                  CGNAMES(ATOMCOUNTER) = "MG"
                  CGTYPE(ATOMCOUNTER) = 16
                  CGMASS(ATOMCOUNTER) = 24.305D0
                  CGCHARGE(ATOMCOUNTER) = 2.000D0
               ELSE IF (CGRESNAMES(J).EQ."NA") THEN
                  CGNAMES(ATOMCOUNTER) = "NA"
                  CGTYPE(ATOMCOUNTER) = 17
                  CGMASS(ATOMCOUNTER) = 22.990D0
                  CGCHARGE(ATOMCOUNTER) = 1.000D0 
               ELSE IF (CGRESNAMES(J).EQ."K") THEN
                  CGNAMES(ATOMCOUNTER) = "K"
                  CGTYPE(ATOMCOUNTER) = 18
                  CGMASS(ATOMCOUNTER) = 39.098D0
                  CGCHARGE(ATOMCOUNTER) = 1.000D0
               ELSE IF (CGRESNAMES(J).EQ."CL") THEN
                  CGNAMES(ATOMCOUNTER) = "CL"
                  CGTYPE(ATOMCOUNTER) = 19
                  CGMASS(ATOMCOUNTER) = 35.450D0
                  CGCHARGE(ATOMCOUNTER) = -1.000D0
               END IF               
            ELSE IF ((RESTYPE(J).EQ.0).OR.(RESTYPE(J).EQ.1)) THEN
               IF (INDEX(CGRESNAMES(J),"5").EQ.0) THEN
                  !first grain: P (not existent in 5 terminal nucleotides)
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "P"
                  CGTYPE(ATOMCOUNTER) = 3
                  CGMASS(ATOMCOUNTER) = 36.97D0
                  CGCHARGE(ATOMCOUNTER) = -1.000D0
               END IF
               !second grain: O5
               ATOMCOUNTER = ATOMCOUNTER + 1
               CGNAMES(ATOMCOUNTER) = "O5"
               CGTYPE(ATOMCOUNTER) = 2
               CGMASS(ATOMCOUNTER) = 16.000D0
               CGCHARGE(ATOMCOUNTER) = 0.000D0
               !third grain: O3
               ATOMCOUNTER = ATOMCOUNTER + 1
               CGNAMES(ATOMCOUNTER) = "O3"
               CGTYPE(ATOMCOUNTER) = 1
               CGMASS(ATOMCOUNTER) = 16.000D0
               CGCHARGE(ATOMCOUNTER) = 0.000D0
               !fourth grain differs between RNA and DNA!
               ATOMCOUNTER = ATOMCOUNTER + 1
               IF (RESTYPE(J).EQ.0) THEN !RNA
                  !third grain: R4
                  CGNAMES(ATOMCOUNTER) = "R4"
                  CGTYPE(ATOMCOUNTER) = 4
               ELSE IF (RESTYPE(J).EQ.1) THEN !DNA
                  !third grain: S4
                  CGNAMES(ATOMCOUNTER) = "S4"
                  CGTYPE(ATOMCOUNTER) = 5
               END IF
               CGMASS(ATOMCOUNTER) = 20.000D0
               CGCHARGE(ATOMCOUNTER) = 0.000D0 
               !fifth grain also differs between RNA and DNA
               ATOMCOUNTER = ATOMCOUNTER + 1
               IF (RESTYPE(J).EQ.0) THEN !RNA
                  !fourth grain: R1
                  CGNAMES(ATOMCOUNTER) = "R1"
                  CGTYPE(ATOMCOUNTER) = 6
               ELSE IF (RESTYPE(J).EQ.1) THEN !DNA
                  !fourth grain: S1
                  CGNAMES(ATOMCOUNTER) = "S1"
                  CGTYPE(ATOMCOUNTER) = 7
               END IF
               CGMASS(ATOMCOUNTER) = 12.000D0
               CGCHARGE(ATOMCOUNTER) = 0.000D0 
               !nucleobase grains
               IF (INDEX(CGRESNAMES(J),"A").GT.0) THEN
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "A1"
                  CGTYPE(ATOMCOUNTER) = 10
                  CGMASS(ATOMCOUNTER) = 67.000D0
                  CGCHARGE(ATOMCOUNTER) = 0.000D0
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "A2"
                  CGTYPE(ATOMCOUNTER) = 11
                  CGMASS(ATOMCOUNTER) = 67.000D0
                  CGCHARGE(ATOMCOUNTER) = 0.000D0
               ELSE IF (INDEX(CGRESNAMES(J),"C").GT.0) THEN
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "C1"
                  CGTYPE(ATOMCOUNTER) = 13
                  CGMASS(ATOMCOUNTER) = 130.000D0
                  CGCHARGE(ATOMCOUNTER) = 0.000D0
               ELSE IF (INDEX(CGRESNAMES(J),"G").GT.0) THEN
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "G1"
                  CGTYPE(ATOMCOUNTER) = 8
                  CGMASS(ATOMCOUNTER) = 75.000D0
                  CGCHARGE(ATOMCOUNTER) = 0.000D0
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "G2"
                  CGTYPE(ATOMCOUNTER) = 9
                  CGMASS(ATOMCOUNTER) = 75.000D0
                  CGCHARGE(ATOMCOUNTER) = 0.000D0
               ELSE IF (INDEX(CGRESNAMES(J),"T").GT.0) THEN
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "T1"
                  CGTYPE(ATOMCOUNTER) = 14
                  CGMASS(ATOMCOUNTER) = 131.000D0
                  CGCHARGE(ATOMCOUNTER) = 0.000D0
               ELSE IF (INDEX(CGRESNAMES(J),"U").GT.0) THEN   
                  ATOMCOUNTER = ATOMCOUNTER + 1
                  CGNAMES(ATOMCOUNTER) = "U1"
                  CGTYPE(ATOMCOUNTER) = 13
                  CGMASS(ATOMCOUNTER) = 131.000D0
                  CGCHARGE(ATOMCOUNTER) = 0.000D0
               END IF
            END IF
         END DO
      END SUBROUTINE ASSIGN_GRAIN_DATA

END MODULE CG_DATA
