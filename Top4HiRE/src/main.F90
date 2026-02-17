PROGRAM CREATE_TOPOLOGY
   USE TOP_GLOBALS
   USE PARSE_AA_PDB, ONLY: PARSE_PDB_INPUT
   !USE PARSE_CG_PDB, ONLY: PARSE_CG_INPUT
   !USE PARSE_SEQ, ONLY: PARSE_FASTA_INPUT
   IMPLICIT NONE
   INTEGER :: NARGS
   INTEGER, PARAMETER :: STDOUT = 6 
   CHARACTER(LEN=50) :: INPUTNAME      !Name of topology file 
   WRITE(STDOUT,'(A)') "Creating Topology file HiRE"
   ! check number of arguments
   NARGS = COMMAND_ARGUMENT_COUNT()
   ! We expect two arguments, the input file and the file type
   IF (NARGS.EQ.2) THEN
      CALL GET_COMMAND_ARGUMENT(1, INPUTNAME)
      CALL GET_COMMAND_ARGUMENT(2, MODE) 
   ELSE
      WRITE(STDOUT,'(A,I4)') "Expecting two arguments, but got ", NARGS
      STOP
   END IF
   ! call input parsing
   IF (MODE.EQ."PDB") THEN
      WRITE(STDOUT,*) " Creating topology from all atom pdb file: ", INPUTNAME
      CALL PARSE_PDB_INPUT(INPUTNAME)
   ELSE IF (MODE.EQ."CG") THEN
      WRITE(STDOUT,*) " Creating topology from CG pdb file: ", INPUTNAME
      !CALL PARSE_CG_INPUT(INPUTNAME)
   ELSE IF (MODE.EQ."SEQ") THEN
      WRITE(STDOUT,*) " Creating topology from fasta file: ", INPUTNAME
      !CALL PARSE_FASTA_INPUT(INPUTNAME)
   ELSE
      WRITE(STDOUT,*) " ERROR - Unknown input mode: ", MODE, " - Needs to be one of CG, PDB or SEQ."
      STOP
   END IF















END PROGRAM CREATE_TOPOLOGY