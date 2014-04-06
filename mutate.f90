SUBROUTINE Find_Mut_Pot_Misprimes()
!
! This is a position dependent replacement for Find_Potential_Misprimes.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,start,finish
  LOGICAL,EXTERNAL :: HMatchNum

  IF (TEST2) PRINT *,"Find_Mut_Pot_Misprimes" !TEST2

! Get rid of potential misprimes in the current mutant range

  CALL Decrement_Misprime_Arrays

! Make sure the search doesn't go beyond the possible ranges

  start=MAX(1,MutNtPos(1)-MPLn)
  finish=MIN((MutNtPos(MutNtNum)+1),(DNAlen-MPLn+1))

! If MutNtPos(1) <= MPLn+1, only run the second half-search

  IF (MutNtPos(1).gt.(MPLn+1)) THEN

! First half-search

    DO i=1,MutNtPos(1)-MPLn-1
      DO j=MutNtPos(1)-MPLn,finish
        IF (HMatchNum(i,j,1)) THEN
        IF (CurrDNA%ntID_Tip(i+MPLn-MPTip).eq.CurrDNA%ntID_Tip(j+MPLn-MPTip)) &
          CALL Increment_Misprime_Arrays(i,j,1)
        IF (CurrDNA%ntID_Tip(i).eq.CurrDNA%ntID_Tip(j)) &
          CALL Increment_Misprime_Arrays(i,j,4)
        END IF
        IF (HMatchNum(i,j,-1)) THEN
        IF (CurrDNA%ntID_Tip(i).eq.CurrDNA%ntID_TipRC(j+MPLn-MPTip)) &
          CALL Increment_Misprime_Arrays(i,j,2)
        IF (CurrDNA%ntID_Tip(i+MPLn-MPTip).eq.CurrDNA%ntID_TipRC(j)) &
          CALL Increment_Misprime_Arrays(i,j,3)
        END IF
      END DO
    END DO
  END IF

! Second half-search

  DO i=start,finish
    DO j=i,DNAlen-MPLn+1
      IF (HMatchNum(i,j,1)) THEN
        IF (CurrDNA%ntID_Tip(i+MPLn-MPTip).eq.CurrDNA%ntID_Tip(j+MPLn-MPTip)) &
          CALL Increment_Misprime_Arrays(i,j,1)
        IF (CurrDNA%ntID_Tip(i).eq.CurrDNA%ntID_Tip(j)) &
          CALL Increment_Misprime_Arrays(i,j,4)
      END IF
      IF (HMatchNum(i,j,-1)) THEN
        IF (CurrDNA%ntID_Tip(i).eq.CurrDNA%ntID_TipRC(j+MPLn-MPTip)) &
          CALL Increment_Misprime_Arrays(i,j,2)
        IF (CurrDNA%ntID_Tip(i+MPLn-MPTip).eq.CurrDNA%ntID_TipRC(j)) &
          CALL Increment_Misprime_Arrays(i,j,3)
      END IF
    END DO
  END DO

  CALL Sort_Misprime_Arrays

END SUBROUTINE Find_Mut_Pot_Misprimes
SUBROUTINE Find_Mutated_Repeats()

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,start,finish
  LOGICAL,EXTERNAL :: PairWithinKnownRepeat

  IF (TEST2) PRINT *,"Find_Mutated_Repeats" !TEST2

  CALL Decrement_Repeat_Arrays

! Make sure the search doesn't go beyond the possible ranges

  start=MAX(1,MutNtPos(1)-RepLen)
  finish=MIN((MutNtPos(MutNtNum)+1),(DNAlen-RepLen+1))

! If MutNtPos(1) <= RepLen+1, only run the second half-search

  IF (MutNtPos(1).gt.(RepLen+1)) THEN

! First half-search

    DO i=1,(MutNtPos(1)-RepLen-1)
      DO j=(MutNtPos(1)-RepLen),finish

! Direct repeat search

        IF (i.ne.j) THEN
          IF (.not.PairWithinKnownRepeat(i,j,1)) THEN
            IF (CurrDNA%ntID_Rep(i).eq.CurrDNA%ntID_Rep(j)) &
              CALL Increment_Repeat_Arrays(i,j,1)
          END IF
        END IF

! Inverse repeat search

        IF (.not.PairWithinKnownRepeat(i,j,-1)) THEN
          IF (CurrDNA%ntID_Rep(i).eq.CurrDNA%ntID_RepRC(j)) &
            CALL Increment_Repeat_Arrays(i,j,-1)
        END IF
      END DO
    END DO
  END IF

! Second half-search

  DO i=start,finish
    DO j=i,DNAlen-RepLen+1

! Direct repeat search

      IF (i.ne.j) THEN
        IF (.not.PairWithinKnownRepeat(i,j,1)) THEN
          IF (CurrDNA%ntID_Rep(i).eq.CurrDNA%ntID_Rep(j)) &
            CALL Increment_Repeat_Arrays(i,j,1)
        END IF
      END IF

! Inverse repeat search

      IF (.not.PairWithinKnownRepeat(i,j,-1)) THEN
        IF (CurrDNA%ntID_Rep(i).eq.CurrDNA%ntID_RepRC(j)) &
          CALL Increment_Repeat_Arrays(i,j,-1)
      END IF

    END DO
  END DO

  CALL Sort_Repeat_Arrays

END SUBROUTINE Find_Mutated_Repeats
SUBROUTINE Mutate_Sequence
!
! Mutate a single codon to an alternate codon.  The residue choice is
! determined in Mutate_Wheel, and the codon choice is made here.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,choice
  INTEGER :: AAN                 ! amino acid choice (1-21)
!  INTEGER :: refCN               ! old codon number (1-64)
  INTEGER :: i1,i2,i3            ! nt positions of choice
  REAL :: rand
!  CHARACTER(LEN=3) :: refCOD     ! old codon sequence

  IF (TEST1) PRINT *,"Mutate_Sequence"

! Generate XScores

  CALL Equalize_Scores

! Determine position to mutate

  CALL Mutate_Wheel

! Determine what amino acid exists at the selected residue MutProtPos.
! AAN is a number between 1 and 21, corresponding to a specific amino acid

  AAN=prot2aa(MutProtPos)

! Determine the nt positions and the actual codon sequence for that residue

  i1=prot2nt(MutProtPos)-1
  i2=prot2nt(MutProtPos)
  i3=prot2nt(MutProtPos)+1
!
! randomly choose codon and make sure it's different and available for that AAN

  IF (AAT(AAN)%NumOfActiveCodons.eq.2) THEN
    choice=1
    IF (AAT(AAN)%Codon(choice).eq.CurrDNA%prot2cod(MutProtPos)) choice=2
  ELSE
    choose: DO i=1,1000
      CALL RANDOM_NUMBER(rand)
      choice=(INT(rand*(AAT(AAN)%NumOfActiveCodons)))+1
      IF (AAT(AAN)%Codon(choice).ne.CurrDNA%prot2cod(MutProtPos)) EXIT choose
    END DO choose
  END IF

 ! update the arrays and change the DNA sequence, then quit

  CurrDNA%prot2cod(MutProtPos)=AAT(AAN)%Codon(choice)
  CurrDNA%nt2cod(i2) = AAT(AAN)%Codon(choice)

  IF (ChainReverse(prot2chain(MutProtPos))) THEN
    CurrDNA%DNAseq(i1:i3)=CFT(AAT(AAN)%Codon(choice))%SeqRC
    CurrDNA%NUMseq(i1)=CFT(AAT(AAN)%Codon(choice))%numRC(1)
    CurrDNA%NUMseq(i2)=CFT(AAT(AAN)%Codon(choice))%numRC(2)
    CurrDNA%NUMseq(i3)=CFT(AAT(AAN)%Codon(choice))%numRC(3)
  ELSE
    CurrDNA%DNAseq(i1:i3)=CFT(AAT(AAN)%Codon(choice))%Seq
    CurrDNA%NUMseq(i1)=CFT(AAT(AAN)%Codon(choice))%num(1)
    CurrDNA%NUMseq(i2)=CFT(AAT(AAN)%Codon(choice))%num(2)
    CurrDNA%NUMseq(i3)=CFT(AAT(AAN)%Codon(choice))%num(3)
  END IF

! Set the MutNtPos values

  MutNtNum=0
  MutNtPos(1)=0
  MutNtPos(2)=0
  MutNtPos(3)=0

  IF (CurrDNA%NUMseq(i1).ne.StoreDNA%NUMseq(i1)) THEN
    MutNtNum=MutNtNum+1
    MutNtPos(MutNtNum)=i1
  END IF
  IF (CurrDNA%NUMseq(i2).ne.StoreDNA%NUMseq(i2)) THEN
    MutNtNum=MutNtNum+1
    MutNtPos(MutNtNum)=i2
  END IF
  IF (CurrDNA%NUMseq(i3).ne.StoreDNA%NUMseq(i3)) THEN
    MutNtNum=MutNtNum+1
    MutNtPos(MutNtNum)=i3
  END IF

  SequenceTranslated=.TRUE.

END SUBROUTINE Mutate_Sequence
SUBROUTINE Mutate_Wheel
!
! Choose a position to mutate based on the score of mutatable codons

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k
  REAL :: rand
  REAL :: ZScore(3333)          ! an accumulated codon-based overall score
  REAL :: choice

  IF (TEST1) PRINT *,"Mutate_Wheel"

! Generate ZScore array

  ZScore(1)=XScore(1)

! If there are more than one codon to be mutated,

  IF (mutPROTnum.gt.1) THEN

    DO i=2,mutPROTnum
      ZScore(i)=ZScore(i-1)+XScore(i)
    END DO

! Pick a random number between 0 and the sum of all the xScore values.
  
    CALL RANDOM_NUMBER(rand)
    choice=rand*ZScore(mutPROTnum)

! Find the codon that corresponds to this number, assign to MutProtPos

    inner: DO j=1,mutPROTnum
      IF (ZScore(j).ge.choice) THEN
        MutProtPos=mutPROT2prot(j)
        EXIT inner
      END IF
    END DO inner
  ELSE

! Otherwise, just choose the first codon

    MutProtPos=mutPROT2prot(1)

  END IF

END SUBROUTINE Mutate_Wheel
