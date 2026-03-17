PROGRAM CREATE_TOPOLOGY
   USE TOP_GLOBALS
   USE FF_GLOBALS, ONLY: DEALLOC_FF_GLOBALS
   USE PARSE_FF, ONLY: PARSE_FF_FILES
   USE PARSE_AA_PDB, ONLY: PARSE_AA_INPUT
   USE PARSE_CG_PDB, ONLY: PARSE_CG_INPUT
   !USE PARSE_SEQ, ONLY: PARSE_FASTA_INPUT
   USE CREATE_TOP, ONLY: WRITE_TOPOLOGY
   USE WRITE_COORDS, ONLY: WRITE_START, WRITE_XYZ
   IMPLICIT NONE
   INTEGER :: NARGS, I
   INTEGER, PARAMETER :: STDOUT = 6
   CHARACTER(LEN=50) :: INPUTNAME      !Name of topology file
   INTEGER, PARAMETER :: NFF = 1
   CHARACTER(LEN=250) :: FFFILES(NFF)
   LOGICAL :: EXISTS

   WRITE(STDOUT,'(A)') "Creating Topology file HiRE"
   ! check number of arguments
   NARGS = COMMAND_ARGUMENT_COUNT()
   ! We expect two arguments, the input file and the file type
   IF (NARGS.EQ.2+NFF) THEN
      CALL GET_COMMAND_ARGUMENT(1, INPUTNAME)
      CALL GET_COMMAND_ARGUMENT(2, MODE)
      DO I=1,NFF
         CALL GET_COMMAND_ARGUMENT(I+2, FFFILES(I))
      END DO
   ELSE
      WRITE(STDOUT,'(A,I4)') "Expecting two+NFF arguments, but got ", NARGS
      STOP
   END IF
   ! call input parsing
   IF (MODE.EQ."PDB") THEN
      WRITE(STDOUT,*) " Creating topology from all atom pdb file: ", INPUTNAME
      CALL PARSE_AA_INPUT(INPUTNAME)
   ELSE IF (MODE.EQ."CG") THEN
      WRITE(STDOUT,*) " Creating topology from CG pdb file: ", INPUTNAME
      CALL PARSE_CG_INPUT(INPUTNAME)
   ELSE IF (MODE.EQ."SEQ") THEN
      WRITE(STDOUT,*) " Creating topology from fasta file: ", INPUTNAME
      !CALL PARSE_FASTA_INPUT(INPUTNAME)
   ELSE
      WRITE(STDOUT,*) " ERROR - Unknown input mode: ", MODE, " - Needs to be one of CG, PDB or SEQ."
      STOP
   END IF
   !parse force field
   CALL PARSE_FF_FILES(NFF,FFFILES)
   !write the topology
   CALL WRITE_TOPOLOGY()
   !write coordinates
   CALL WRITE_START()
   CALL WRITE_XYZ()

   !tidy up
   CALL DEALLOC_FF_GLOBALS()
   CALL DEALLOC_TOP_GLOBALS()
END PROGRAM CREATE_TOPOLOGY
