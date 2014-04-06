SUBROUTINE Fix_Degenerates
!
! Fix degenerate sequences to A,C,G or T

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k
  REAL :: rand

  IF (TEST0) PRINT *,"Fix_Degenerates" !TEST0

  DO i=1,NumDegPos

! Choose a nt at random from possible
    CALL RANDOM_NUMBER(rand)
    j = CurrDNA%DegenNum(DegPos(i)) ! index for degenerate sequence
    k=(INT(rand*(DegenSeq(j)%NumOfNT)))+1 ! choice for index
    CurrDNA%DNAseq(DegPos(i):DegPos(i)) = DegenSeq(j)%Seq(k) ! assign seq
    CurrDNA%NumSeq(DegPos(i)) = DegenSeq(j)%NumSeq(k) ! assign NumSeq
  
  END DO

END SUBROUTINE Fix_Degenerates
SUBROUTINE IntSwap(firstelement,lastelement)

! Swap integers

  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: firstelement,lastelement,dummy

  IF (TEST3) PRINT *,"IntSwap" !TEST3

  dummy=firstelement
  firstelement=lastelement
  lastelement=dummy

END SUBROUTINE IntSwap
SUBROUTINE RealSwap(firstelement,lastelement)

! Swap real numbers

  USE dnaworks_test
  IMPLICIT NONE

  REAL :: firstelement,lastelement,dummy

  IF (TEST3) PRINT *,"RealSwap" !TEST3

  dummy=firstelement
  firstelement=lastelement
  lastelement=dummy

END SUBROUTINE RealSwap
SUBROUTINE Revert_Degenerates
!
! Revert degenerate sequences back to original

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j

  IF (TEST0) PRINT *,"Revert_Degenerates" !TEST0

  DO i=1,NumDegPos

    j = CurrDNA%DegenNum(DegPos(i)) ! index for degenerate sequence at that position
    CurrDNA%DNAseq(DegPos(i):DegPos(i)) = DegenSeq(j)%DegNT ! assign seq
  
  END DO

END SUBROUTINE Revert_Degenerates
