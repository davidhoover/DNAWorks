PROGRAM dnaworks

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER,EXTERNAL :: CurrentTimeSeconds

  IF (TEST0) PRINT *,"DNAWORKS start" !TEST0
  CALL RANDOM_SEED()

  MainTimeStart=CurrentTimeSeconds()  ! when does the run begin?

  CALL Get_Args                 ! get the command arguments, if any
  CALL Default_Param
  CALL Read_Input

  IF (DNAlen.LE.50) CALL Stop_Program("DNA length is less than 50 nt.")

! Reset weights if there is no protein to mutate

  IF (PROTlen.eq.0) THEN
    Cwt=0.0                       ! weight for codon scoring
    Rwt=0.0                       ! weight for repeat scoring
    Gwt=0.0                       ! weight for GC scoring
    Awt=0.0                       ! weight for AT scoring
    Pwt=0.0                       ! weight for pattern scoring
  END IF

! start logfile

  OPEN (UNIT=outputnum,FILE=outputfile,FORM="FORMATTED",STATUS="REPLACE")

  CALL Print_Output_Start(outputnum)
  CALL Print_Output_Start(console)
  CALL FLUSH(console)
  CALL Print_Seq_Log(outputnum)
  IF (PROTlen.gt.0) CALL Print_Codon_Log(outputnum)
  CALL Print_Pattern_Log(outputnum)
!  CALL Print_TranslatedDNA(outputnum)

! determine the number of solutions

  TotalNumberOfSolutions=(NumberOfSolutions*(MeltTempHi-MeltTempLo+1)* &
      & (OligoLenHi-OligoLenLo+1))
  IF (TotalNumberOfSolutions.gt.9999) CALL Stop_Program("Too many trials.  Limit the range of parameters.")
  WRITE(UNIT=console,FMT="('')")
  WRITE(UNIT=console,FMT="(20x,'Starting ',i3,' trial',$)") TotalNumberOfSolutions
  IF (TotalNumberOfSolutions.gt.1) WRITE(UNIT=console,FMT="('s',$)")
  WRITE(UNIT=console,FMT="('...')")
  CALL FLUSH(console)

  CALL Run_Dnaworks()

  CALL Print_Final_Tally(console)
  CALL Print_Output_End(console)
  CALL FLUSH(console)
  CALL Print_Final_Tally(outputnum)
  CALL Print_Output_End(outputnum)

  CLOSE (UNIT=outputnum)

  IF (LEN_TRIM(email).GT.5) THEN
    CALL Send_Email
  END IF

  WRITE(UNIT=console,FMT="(' ')")
  WRITE(UNIT=console,FMT="('Finished ')")

END PROGRAM dnaworks
