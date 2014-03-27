SUBROUTINE Get_Args

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=30) :: ARGV(100) ! command line arguments
  INTEGER :: ARGC                 ! number of command line arguments
  INTEGER :: i,j

! IARGC returns the total number of arguments on the command line

  ARGC=IARGC()

! GETARG returns the argument that corresponds to the argument number, with
! zero equal to the command itself

  DO i=1,ARGC
    CALL GETARG(i,ARGV(i))

! Turn on testing mode

    IF (INDEX(ARGV(i),"-t3").eq.1) THEN
      TEST3=.TRUE.
      TEST2=.TRUE.
      TEST1=.TRUE.
      TEST0=.TRUE.
    ELSE IF (INDEX(ARGV(i),"-t2").eq.1) THEN
      TEST2=.TRUE.
      TEST1=.TRUE.
      TEST0=.TRUE.
    ELSE IF (INDEX(ARGV(i),"-t1").eq.1) THEN
      TEST1=.TRUE.
      TEST0=.TRUE.
    ELSE IF (INDEX(ARGV(i),"-t0").eq.1) THEN
      TEST0=.TRUE.
    ELSE IF (INDEX(ARGV(i),"-q").eq.1) THEN
      QUIET=.TRUE.
    ELSE IF (INDEX(ARGV(i),"-fast").eq.1) THEN
      FAST=.TRUE.
    ELSE IF (INDEX(ARGV(i),"-help").eq.1) THEN
      CALL Print_Help
    ELSE

! Assign inputfile from ARGV(1) if possible

      inputfile=ARGV(i)
    END IF
  END DO

END SUBROUTINE Get_Args
SUBROUTINE Oligo_Design(SolutionNo,num)
!
! This subroutine is the main engine of the program.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,main_count,z,timediff
  INTEGER :: nlimit  ! max number of successful changes before continuing
  INTEGER :: nover  ! max number of changes before dropping the temperature
  INTEGER :: maxlimit  ! max number of counts before quitting
  INTEGER :: nsucc
  INTEGER :: count
  INTEGER :: SolutionNo ! current solution number
  INTEGER :: num        ! where to print output?
  LOGICAL :: ans
  REAL :: t         ! initial temperature
  REAL :: tfactr   ! how much to drop temperature
  REAL :: gain
  REAL :: guess
  REAL :: rand
  INTEGER,EXTERNAL :: CurrentTimeSeconds

  IF (TEST0) PRINT *,"Oligo_Design" !TEST0

  nlimit=50
  nover=500
  nsucc=0
  count=0
  t=0.5
  tfactr=0.96
  main_count=0
  maxlimit=1500
  IF (FAST) maxlimit=300

  WRITE(UNIT=console,FMT="('')")

  BestDNA = CurrDNA     ! initialize BestDNA

  tempdrop: DO i=1,1000 ! There will be a total of 1000 drops in temperature

! just dump out if there are no protein residues

    IF (PROTlen.eq.0.and.(.not.OligoLenRandom)) THEN
      BestDNA = CurrDNA
      EXIT tempdrop
    END IF

    nsucc=0

    mutate: DO j=1,nover ! within which there will be <nover> mutation/length change rounds.

      StoreDNA = CurrDNA ! Make a backup of the current solution

      IF (PROTlen.gt.0) CALL Mutate_Sequence  ! only mutate protein sequence
      main_count=main_count+1

      CALL Generate_Overlaps(SolutionNo) ! Generate overlaps for the mutated sequence
  
      IF (MOD(CurrDNA%NumOlaps,2).eq.0) CALL Stop_Program("Even number of overlaps.  Try adjusting parameters.")

      gain = CurrDNA%OverallScore - StoreDNA%OverallScore ! Determine whether the mutated solution is any better than the old one

! If the gain is good enough or if the temperature is high enough, we have a successful mutation round.

      CALL RANDOM_NUMBER(rand) 
      ans=(gain.lt.0.0).or.(rand.lt.exp(-gain/t))        !  Metropolis
      IF (ans) THEN
        nsucc=nsucc+1
      ELSE
        CurrDNA=StoreDNA ! If not, go back to the original sequence and try again.
      END IF

      IF (CurrDNA%OverallScore.lt.BestDNA%OverallScore) THEN
        BestDNA = CurrDNA ! If the current sequence is better than the best sequence, replace the BestDNA with CurrDNA.  
        count=0
      ELSE
        count=count+1 ! If the current score does not achieve a better value than the best score, then start counting.
      END IF

      IF ((MOD(main_count,100)).eq.0) WRITE(UNIT=num,FMT=&
        "(6x,i5,' optimization rounds, best = ',f9.3,' Rep =',i4,' Mis =',i4)") &
        main_count,BestDNA%OverallScore,BestDNA%RN,BestDNA%MSN ! Keep the user informed

      IF (count.gt.maxlimit) THEN ! If the count between drops is greater than 300, then quit.
        EXIT tempdrop
      END IF

      IF (nsucc.ge.nlimit) EXIT mutate ! If there are more than <nlimit> successful rounds, exit the mutation loop and move to the next temperature drop.

      timediff = CurrentTimeSeconds()-MainTimeStart;
      IF (MainTimeLimit.GT.0.and.timediff.GE.MainTimeLimit) THEN ! Dump out if out of time
        WRITE(UNIT=console,FMT="(/,'Main time limit reached.')")
        WRITE(UNIT=outputnum,FMT="(/,'Main time limit reached.')")
        TimesUp=.TRUE.
        EXIT mutate
      END IF

      CALL FLUSH(console)

    END DO mutate

    t=t*tfactr ! Drop the temperature

    IF (nsucc.eq.0.or.t.lt.0.0001.or.TimesUp) THEN ! If within <nover> rounds of mutation there are no successes or the temperature is too low, or out of time, then quit.  
      WRITE(UNIT=num,FMT="('Limit of simulated annealing, quitting.')")
      EXIT tempdrop
    END IF

    IF (BestDNA%OverallScore.lt.0.001) EXIT tempdrop ! If the best score is good enough, then quit.

  END DO tempdrop

  CurrDNA = BestDNA ! Push the best solution from Oligo_Design into the current solution

  CALL Revert_Degenerates
  CALL Print_FinalDNA_Log(outputnum,SolutionNo)
  CALL Print_Scores_Log(console)
  CALL Print_Scores_Log(outputnum)
  CALL Print_Histogram(outputnum,SolutionNo)
  CALL Print_Pattern_Screen(outputnum)
  CALL Print_Oligo_Log(outputnum)

END SUBROUTINE Oligo_Design
SUBROUTINE Run_Dnaworks()

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,l,timediff
  INTEGER,EXTERNAL :: CurrentTimeSeconds

  IF (TEST0) PRINT *,"Run_Dnaworks" !TEST0

! start the loops

  melt: DO j=MeltTempLo,MeltTempHi
    MeltTemp=j
    oligo: DO l=OligoLenLo,OligoLenHi
      OligoLen=l
      main: DO i=1,NumberOfSolutions
        SequenceTranslated=.FALSE.
        CALL Translate_Protein

! Dump out if out of time

        timediff = CurrentTimeSeconds()-MainTimeStart;
        IF (MainTimeLimit.GT.0.and.timediff.GE.MainTimeLimit) THEN
          WRITE(UNIT=console,FMT="(/,'Main time limit reached.')")
          WRITE(UNIT=outputnum,FMT="(/,'Main time limit reached.')")
          TimesUp=.TRUE.
          FinalScore(CurrSolutionNo)%Oligo=OligoLen
          FinalScore(CurrSolutionNo)%MeltT=MeltTemp
          EXIT melt
        END IF

        CurrSolutionNo = CurrSolutionNo+1

        CALL Print_Param_Log(console,CurrSolutionNo)
        CALL Print_Param_Log(outputnum,CurrSolutionNo)
        CALL FLUSH(console)

        CALL Generate_Overlaps(CurrSolutionNo)
  
        IF (MOD(CurrDNA%NumOlaps,2).eq.0) THEN
          IF (.not.QUIET) THEN
            WRITE(UNIT=console,FMT="('Even number of overlaps - trial ',i4,' abandoned')") CurrSolutionNo
            WRITE(UNIT=outputnum,FMT="('Even number of overlaps - trial ',i4,' abandoned')") CurrSolutionNo
          END IF
          FinalScore(CurrSolutionNo)%FinaScore=999999
          CYCLE main
        END IF

! If everything is ok, go to Oligo_Design

        CALL Oligo_Design(CurrSolutionNo,console)

! Keep track of times

        CALL Print_Estimated_Time(CurrSolutionNo)

! Update FinalScore tally

        FinalScore(CurrSolutionNo)%Oligo=OligoLen
        FinalScore(CurrSolutionNo)%MeltT=MeltTemp

      END DO main
    END DO oligo
  END DO melt

! in case optimization stopped prematurely

  TotalNumberOfSolutions=CurrSolutionNo

END SUBROUTINE Run_Dnaworks
SUBROUTINE Stop_Program(message)
!
! This subroutine stops the program and displays an error message.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: message
  INTEGER :: i

  IF (TEST0) PRINT *,"Stop_Program" !TEST0

  WRITE(UNIT=console,FMT="(' ')")
  WRITE(UNIT=console,FMT="('Program error:')")
  WRITE(UNIT=console,FMT="(a)") message
  WRITE(UNIT=console,FMT="('Exiting program now')")
  CALL FLUSH(console)

  WRITE(UNIT=outputnum,FMT="(' ')")
  WRITE(UNIT=outputnum,FMT="('Program error:')")
  WRITE(UNIT=outputnum,FMT="(a)") message
  WRITE(UNIT=outputnum,FMT="('Exiting program now')")

  CLOSE (UNIT=outputnum)
  STOP

END SUBROUTINE Stop_Program
