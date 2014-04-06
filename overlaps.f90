INTEGER FUNCTION ForOlap(first)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: first,last     ! positions in DNAseq
  REAL,EXTERNAL :: TmCalc
  REAL :: diff,diff_lo,diff_hi,diff2
  LOGICAL :: done
  INTEGER :: shift

  IF (TEST3) PRINT *,"ForOlap" !TEST3

  done=.FALSE.
  shift=32

  last=first+shift
  shift=shift/2

  loop: DO WHILE (.not.done)
    IF (last.ge.DNAlen) THEN
      last=last-shift
    ELSE
      diff=MeltTemp-(TmCalc(first,last))
      IF (ABS(diff).gt.MeltTol) THEN
        IF (diff.gt.0) THEN
          last=last+shift
        ELSE
          last=last-shift
        END IF
      ELSE
        done=.TRUE.
      END IF
    END IF
    shift=shift/2
    IF (shift.le.1) EXIT loop
  END DO loop

! For the final step, determine which of the final two positions is best

  IF (.not.done) THEN
    IF (last.le.(DNAlen-1)) THEN
      diff=MeltTemp-(TmCalc(first,last))
      IF (diff.gt.0) THEN
        shift=1
      ELSE
        shift=-1
      END IF
      last=last+shift
      diff2=MeltTemp-(TmCalc(first,last))
      IF (ABS(diff).lt.ABS(diff2)) last=last-shift
    END IF
    IF ((DNAlen-last).le.2) last=DNAlen
  END IF

  ForOlap = last

END FUNCTION ForOlap
SUBROUTINE Generate_Overlaps(SolutionNo)
!
! The nucleotide sequence is broken into overlaps of around 20 nucleotides
!   each, depending on the calculated Tm.  The set of potential oligos is
!   then analyzed and the best trial is kept.  A gap is allowed between
!   overlaps to give oligos of size oligoLen.
!
! The structure of the overlap array is as follows
!
!         1,1   1,2       2,1     2,2   3,1  3,2        4,1  4,2   5,1   5,2
! ................         ....................          ..................
!          .........................     ......................     ..........
!          -------         ---------     ------          ------     -------
!
!   OVERLAP:  1                2           3                4          5
!

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,l,m
  INTEGER :: first          ! first nt of the overlap
  INTEGER :: last=1         ! last nt of the overlap
  INTEGER :: SolutionNo     ! current solution number
  INTEGER :: shift          ! number of nt before first overlap
  INTEGER :: reset          ! first overlap size
  REAL :: rand
  INTEGER :: olength        ! the number of nt to skip ahead
  LOGICAL :: changed        ! true if CurrDNA%OverallScore.lt.BestOverlapDNA%OverallScore

  IF (TEST1) PRINT *,"Generate_Overlaps"

  changed=.FALSE.

  CALL Fix_Degenerates ! pin down degenerate sequences

  IF (.not.MutantRun) THEN

    generate: DO k=1,10000
  
! initialize the nt2overlap array

      shift=0
      last=1
      changed=.FALSE.

      DO i=1,DNAlen
        nt2overlap=0
      END DO
  
      BestOverlapDNA = CurrDNA            ! initialize BestOverlapDNA values
      BestOverlapDNA%OverallScore = 9999

      DO i=1,999                        ! initialize the arrays
        CurrDNA%OlapsPos(i,1)=0
        CurrDNA%OlapsPos(i,2)=0
      END DO
      CurrDNA%NumOlaps=0
    
! Determine the new CurrDNA%OlapsPos values
  
      outer: DO i=1,1000           ! keep shifting
  
        IF (NOGAPS) THEN
          olength = 0
        ELSE 
          IF (OligoLenRandom) THEN
            CALL RANDOM_NUMBER(rand) 
            rand = (rand*(OligoLen-20))
            olength = (INT(rand))+20        ! randomize oligo length
          ELSE
            olength = OligoLen
          END IF
        END IF
  
        CurrDNA%NumOlaps=0                 ! initialize the number of overlaps
        first=1+shift
        CALL Make_Olap(first,last)
        IF ((shift.gt.0).and.(CurrDNA%OlapsPos(1,2).ge.OligoLen)) EXIT outer
        last=first+olength-1
        first=last-7                           ! the minimal overlap size is 7
        IF (first.le.CurrDNA%OlapsPos(1,2)) THEN
          first = CurrDNA%OlapsPos(1,2)+1
          last = first+7
        END IF
    
        inner: DO j=1,999
  
          IF (NOGAPS) THEN
            olength = 0
          ELSE 
            IF (OligoLenRandom) THEN
              CALL RANDOM_NUMBER(rand) 
              rand = (rand*(OligoLen-20))
              olength = (INT(rand))+20       ! randomize oligo length
            ELSE
              olength = OligoLen
            END IF
          END IF
  
          IF (last.ge.DNAlen) EXIT inner
          CALL Make_Olap(first,last)
          last=first+olength-1
          first=last-7
          IF (first.le.CurrDNA%OlapsPos(CurrDNA%NumOlaps,2)) THEN
            first = CurrDNA%OlapsPos(CurrDNA%NumOlaps,2)+1
            last = first+7
          END IF
        END DO inner
    
        shift=shift+1
    
        IF ((MOD(CurrDNA%NumOlaps,2)).eq.0) THEN
          CYCLE outer
        END IF
    
        CALL Evaluate_Scores
    
        IF (CurrDNA%OverallScore.lt.BestOverlapDNA%OverallScore) THEN
          BestOverlapDNA = CurrDNA
          changed=.TRUE.
        END IF
    
      END DO outer
  
      CurrDNA=BestOverlapDNA                ! revert to best solution
  
      IF (MOD(CurrDNA%NumOlaps,2).eq.1) THEN
        EXIT generate
      ELSE
        IF (TEST0) PRINT *,k,"EVEN OVERLAPS" !TEST0

! Take drastic action to get optimization moving

        IF ((MOD(k,200)).eq.0) THEN
          OligoLenLo=OligoLenLo+1
          OligoLenHi=OligoLenHi+1
          OligoLen=OligoLen+1
          WRITE(UNIT=console,FMT="('')")
          WRITE(UNIT=outputnum,FMT="('')")
          WRITE(UNIT=console,FMT="(' Too many sets of even overlaps -- increasing oligo length to',i4)") OligoLen
          WRITE(UNIT=outputnum,FMT="(' Too many sets of even overlaps -- increasing oligo length to',i4)") OligoLen
        END IF
      END IF
    END DO generate

  END IF

! Assign the nt2overlap array
  
  DO i=1,DNAlen
    DO j=1,CurrDNA%NumOlaps
      IF (i.ge.CurrDNA%OlapsPos(j,1).and.i.le.CurrDNA%OlapsPos(j,2)) THEN
        nt2overlap(i)=j
      END IF
    END DO
  END DO

  IF (.not.changed) CALL Evaluate_Scores ! in case the CurrDNA is never better than BestOverlapDNA

  FinalScore(SolutionNo)%InitScore=CurrDNA%OverallScore
  
END SUBROUTINE Generate_Overlaps
SUBROUTINE Make_Olap(first,last)
!
! Simplifies the process of finding overlaps.  The OlapsPos and MeltT values
! are recorded in this subroutine for each overlap.  It also automates the
! decision making about forward or reverse methods of generating overlaps.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: first,last
  INTEGER,EXTERNAL :: ForOlap
  INTEGER,EXTERNAL :: RevOlap
  REAL,EXTERNAL :: TmCalc

  IF (TEST2) PRINT *,"Make_Olap" !TEST2

  IF (CurrDNA%NumOlaps.eq.0) THEN

    last=ForOlap(first)

    CurrDNA%NumOlaps=1
    CurrDNA%OlapsPos(1,1)=first
    CurrDNA%OlapsPos(1,2)=last
    CurrDNA%MeltT(1)=TmCalc(first,last)

!  PRINT *,first,last,CurrDNA%MeltT(1)
  ELSE

    first=RevOlap(last)

    IF (first.le.CurrDNA%OlapsPos(CurrDNA%NumOlaps,2)) THEN
      first = CurrDNA%OlapsPos(CurrDNA%NumOlaps,2)+1
      last = ForOlap(first)
    END IF

    IF (last.lt.DNAlen) THEN
      CurrDNA%NumOlaps=CurrDNA%NumOlaps+1
      CurrDNA%OlapsPos(CurrDNA%NumOlaps,1)=first
      CurrDNA%OlapsPos(CurrDNA%NumOlaps,2)=last
      CurrDNA%MeltT(CurrDNA%NumOlaps)=TmCalc(first,last)
    END IF
  END IF

END SUBROUTINE Make_Olap
INTEGER FUNCTION RevOlap(last)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: first,last     ! positions in DNAseq
  REAL,EXTERNAL :: TmCalc
  REAL :: diff,diff_lo,diff_hi,diff2
  LOGICAL :: done
  INTEGER :: shift

  IF (TEST3) PRINT *,"RevOlap" !TEST3

  shift=32
  done=.FALSE.

  first=last-shift
  shift=shift/2

  loop: DO WHILE (.not.done)
    IF (first.le.1) THEN
      first=first+shift
    ELSE
      diff=MeltTemp-(TmCalc(first,last))
      IF (ABS(diff).gt.MeltTol) THEN
        IF (diff.gt.0) THEN
          first=first-shift
        ELSE
          first=first+shift
        END IF
      ELSE
        done=.TRUE.
      END IF
    END IF
    shift=shift/2
    IF (shift.le.1) EXIT loop
  END DO loop

! For the final step, determine which of the final two positions is best

  IF (.not.done) THEN
    IF (first.ge.2) THEN
      diff=MeltTemp-(TmCalc(first,last))
      IF (diff.gt.0) THEN
        shift=-1
      ELSE
        shift=1
      END IF
      first=first+shift
      diff2=MeltTemp-(TmCalc(first,last))
      IF (ABS(diff).lt.ABS(diff2)) first=first-shift
    END IF
    IF (first.le.2) first=1
  END IF

  RevOlap = first

END FUNCTION RevOlap
