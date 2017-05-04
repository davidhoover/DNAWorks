SUBROUTINE Create_Final_Arrays
!
! Fill in nt2Solig and nt2Aolig arrays.  Note that we are now using
! a different numbering system than previous version.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i
  INTEGER :: olap           ! overlap number
  LOGICAL :: last           ! last overhang

  IF (TEST0) PRINT *,"Create_Final_Arrays" !TEST0

! Assign the nt2Solig and nt2Aolig arrays

  olap=1
  last=.FALSE.

  DO i=1,MaxDNAlen
    nt2Solig(i)=0
    nt2Aolig(i)=0
  END DO

  DO i=1,DNAlen

! If i goes past the end of an overlap, increment the olap number

    IF (i.eq.(CurrDNA%OlapsPos(olap,2)+1)) olap=olap+1

! Don't let the olap number go above NumOlaps number

    IF (olap.gt.CurrDNA%NumOlaps) THEN
      olap=olap-1
      last=.TRUE.
    END IF

! Insert numbers into the array.  There are four places to be within the
! synthetic gene:

! Downstream of an odd overlap

    IF ((MOD(olap,2).eq.1).and.(i.lt.CurrDNA%OlapsPos(olap,1))) THEN
      nt2Solig(i)=olap
      nt2Aolig(i)=0
    END IF

! Downstream of an even overlap (last overhang is special case)

    IF (last) THEN
      nt2Solig(i)=0
      nt2Aolig(i)=olap+1
    ELSE
      IF ((MOD(olap,2).eq.0).and.(i.lt.CurrDNA%OlapsPos(olap,1))) THEN
        nt2Solig(i)=0
        nt2Aolig(i)=olap
      END IF
    END IF

! In the middle of an odd overlap

    IF ((MOD(olap,2).eq.1).and.(i.ge.CurrDNA%OlapsPos(olap,1)).and.(i.le.CurrDNA%OlapsPos(olap,2))) THEN
      nt2Solig(i)=olap
      nt2Aolig(i)=olap+1
    END IF

! In the middle of an even overlap

    IF ((MOD(olap,2).eq.0).and.(i.ge.CurrDNA%OlapsPos(olap,1)).and.(i.le.CurrDNA%OlapsPos(olap,2))) THEN
      nt2Solig(i)=olap+1
      nt2Aolig(i)=olap
    END IF

  END DO

END SUBROUTINE Create_Final_Arrays
SUBROUTINE Export_Misprime_Pairs(num)

! 1. direct-sense(DS): forward primer mispriming on the sense strand
!
!      -------------->               -------------->
!      |||||||||||||||               ..........|||||
! -------------------------------------------------------
!
! 2. inverse-sense(IS): reverse primer mispriming on the sense strand
! NOTE THAT IF THE FORWARD OLIGO MATCHES M2, MSX = 5, NOT 2
!
!                                    -------------->
!                                    ..........|||||
! -------------------------------------------------------
!      |||||||||||||||
!      <--------------
!
! 3. inverse-antisense(IA): forward primer mispriming on the antisense strand
! NOTE THAT IF THE REVERSE OLIGO MATCHES M2, MSX = 6, NOT 3
!
!      -------------->
!      |||||||||||||||
! -------------------------------------------------------
!                                    |||||..........
!                                    <--------------
!
! 4. direct-antisense(DA): reverse primer mispriming on the antisense strand
!
! -------------------------------------------------------
!      |||||||||||||||               |||||..........
!      <--------------               <--------------
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: num,i,j,k,a
  CHARACTER(LEN=9999) :: text1
  CHARACTER(LEN=9999) :: text2
  CHARACTER(LEN=9999) :: text3
  INTEGER :: ms1,ms2,me1,me2,fract,oligo_num,midpoint
  CHARACTER(LEN=2) :: id(6)

  IF (TEST0) PRINT *,"Export_Misprime_Pairs" !TEST0

  id(1)='DS' ! always forward
  id(2)='IS' ! reverse
  id(3)='IA' ! forward
  id(4)='DA' ! always reverse
  id(5)='IS' ! forward
  id(6)='IA' ! reverse

  midpoint=(CurrDNA%NumOlaps+1)/2

  IF (CurrDNA%MSN.gt.0) THEN

    WRITE(UNIT=num,FMT="(' ')")
    WRITE(UNIT=num,FMT="('  There are ',i3,' potential misprimings with <= ',i2,' non-identical nts:')") CurrDNA%MSN,MaxNonId

    WRITE(UNIT=num,FMT='(1x,a64)') bar64

    text1="  Oligo Type 5'-start     Sequence     3'-start  Identical"
    WRITE(UNIT=num,FMT="(a70)") text1
    WRITE(UNIT=num,FMT='(1x,a64)') bar64
    DO i=1,CurrDNA%MSN
      ms1=CurrDNA%MS1(i)
      ms2=CurrDNA%MS2(i)
      me1=CurrDNA%MS1(i)+MPLn-1
      me2=CurrDNA%MS2(i)+MPLn-1
      text1(1:MPLn) = CurrDNA%DNAseq(ms1:me1)
      text3(1:MPLn) = CurrDNA%DNAseq(ms2:me2)
      fract=0
      oligo_num=0

      IF (MOD(CurrDNA%MSX(i),2).eq.0) THEN
        text1(1:MPLn) = CurrDNA%DNAseq(ms1:me1)
        CALL RevComplStr(text1(1:MPLn))
      END IF

      IF (CurrDNA%MSX(i).eq.3.or.CurrDNA%MSX(i).eq.4.or.CurrDNA%MSX(i).eq.5) THEN
        text3(1:MPLn) = CurrDNA%DNAseq(ms2:me2)
        CALL RevComplStr(text3(1:MPLn))
      END IF

! show the mismatches by lowercase

      DO j=1,MPLn
        IF (text1(j:j).ne.text3(j:j)) THEN
          CALL ToLowerCase(text3(j:j))
          text2(j:j)=' '
        ELSE
          text2(j:j)='|'
          fract=fract+1
        END IF
      END DO

! Figure out which oligo the misprime belongs to

      IF (MOD(CurrDNA%MSX(i),2).eq.1) THEN
        oligo_num=nt2Solig(ms1+MPLn-1)
      ELSE
        oligo_num=nt2Aolig(ms1)
      END IF

      IF (TBIO) THEN
        IF (CurrDNA%MOL(i).lt.midpoint) THEN
          oligo_num=CurrDNA%MOL(i)
        ELSE IF (CurrDNA%MOL(i).gt.midpoint) THEN
          oligo_num=CurrDNA%MOL(i)+1
        END IF
      END IF

      IF (MOD(CurrDNA%MSX(i),2).eq.1) THEN
        WRITE(UNIT=num,FMT="(3x,i3,3x,a2,3x,i4,2x,a,2x,i4,6x,i2,'/',i2)") &
          oligo_num,id(CurrDNA%MSX(i)),ms1,text1(1:MPLn),me1,fract,MPLn
      ELSE
        WRITE(UNIT=num,FMT="(3x,i3,3x,a2,3x,i4,2x,a,2x,i4,6x,i2,'/',i2)") &
          oligo_num,id(CurrDNA%MSX(i)),me1,text1(1:MPLn),ms1,fract,MPLn
      END IF

! show the identical matches

      WRITE(UNIT=num,FMT="(20x,a)") text2(1:MPLn)

! print the other position with numbers

      IF (CurrDNA%MSX(i).eq.1.or.CurrDNA%MSX(i).eq.2.or.CurrDNA%MSX(i).eq.6) THEN
        WRITE(UNIT=num,FMT="(14x,i4,2x,a,2x,i4)") ms2,text3(1:MPLn),me2
      ELSE
        WRITE(UNIT=num,FMT="(14x,i4,2x,a,2x,i4)") me2,text3(1:MPLn),ms2
      END IF

      WRITE(UNIT=num,FMT="('')")
    END DO
    WRITE(UNIT=num,FMT='(1x,a64)') bar64
  END IF

END SUBROUTINE Export_Misprime_Pairs
SUBROUTINE Export_Repeat_Pairs(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: num,i,j,k,a
  CHARACTER(LEN=2) :: txt ! DR, IR, PR
  CHARACTER(LEN=9999) :: RCtext
  INTEGER :: s1,s2,ln,dr
  LOGICAL :: degenWarning

  IF (TEST0) PRINT *,"Export_Repeat_Pairs" !TEST0

  IF (CurrDNA%RN.gt.0) THEN

    WRITE(UNIT=num,FMT="(' ')")
    WRITE(UNIT=num,FMT="('  There are ',i3,' repeats greater than '&
      &,i2,' nt:')")CurrDNA%RN,RepLen

    WRITE(UNIT=num,FMT='(1x,a64)') bar64
    DO i=1,CurrDNA%RN

! Assign the temporary values

      degenWarning = .FALSE.
      s1=CurrDNA%RS1(i)
      s2=CurrDNA%RS2(i)
      ln=CurrDNA%RLn(i)
      dr=CurrDNA%RX(i)

      IF (dr.eq.1) THEN
        txt="DR"
        IF (CurrDNA%DNAseq(s1:(s1+ln-1)).ne.CurrDNA%DNAseq(s2:(s2+ln-1))) degenWarning = .TRUE. 
      ELSE
! Convert second position to RC positions

        RCtext(1:ln) = CurrDNA%DNAseq(s2:(s2+ln-1))
        CALL RevComplStr(RCtext(1:ln))
        IF (CurrDNA%DNAseq(s1:(s1+ln-1)).ne.RCtext(1:ln)) degenWarning = .TRUE.
!          DO j=1,ln
!            CALL ToLowerCase(RCtext(j:j))
!          END DO

        IF (s1.ne.s2) THEN
          txt="IR"
        ELSE
          txt="PR"
        END IF
      END IF
      IF (degenWarning) THEN
        WRITE(UNIT=num,FMT="(1x,a2,1x,'Pos1 = ',i4,3x,'Pos2 = ',i4,3x,'Size = ',i2,3x,'Seq1 = ',a,3x,'WARNING!')") &
        txt(1:2),s1,s2,ln,CurrDNA%DNAseq(s1:(s1+ln-1))
        WRITE(UNIT=num,FMT="(44x,'Seq2 = ',a,3x,'WARNING!')")  CurrDNA%DNAseq(s2:(s2+ln-1))
      ELSE
        WRITE(UNIT=num,FMT="(1x,a2,1x,'Pos1 = ',i4,3x,'Pos2 = ',i4,3x,'Size = ',i2,3x,'Seq1 = ',a)") &
        txt(1:2),s1,s2,ln,CurrDNA%DNAseq(s1:(s1+ln-1))
        WRITE(UNIT=num,FMT="(44x,'Seq2 = ',a)")  CurrDNA%DNAseq(s2:(s2+ln-1))
      END IF
    END DO
    WRITE(UNIT=num,FMT='(1x,a64)') bar64
  END IF

END SUBROUTINE Export_Repeat_Pairs
SUBROUTINE Print_Codon_Histogram(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,num,itmp,CF(12),ctmp(64)
  REAL :: tmp

  IF (TEST0) PRINT *,"Print_Codon_Histogram" !TEST0

! For the logfile: only if protein present


! Detailed Codon Frequency Log

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="('                  DETAILED CODON FREQUENCY REPORT ')")
  WRITE(UNIT=num,FMT="('  [ Codon, AA, Frequency, # of times used in coding sequence ] ')")
  WRITE(UNIT=num,FMT="(1x,a64)") bar64
  WRITE(UNIT=num,FMT="(' ')")

! Detailed Codon use histogram:

  DO i=1,64
    ctmp(i)=0
  END DO

  DO i=1,PROTlen
    ctmp(CurrDNA%prot2cod(i))=ctmp(CurrDNA%prot2cod(i))+1
  END DO

  DO i=1,16
    j=((i-1)*4)+1
    WRITE(UNIT=num,FMT="(1x,a3,1x,a1,1x,f4.2,1x,i3,2x,a3,1x,a1,1x,f4.2,&
      &1x,i3,2x,a3,1x,a1,1x,f4.2,1x,i3,2x,a3,1x,a1,1x,f4.2,1x,i3)") &
       CFT(j)%Seq,CFT(j)%AA1,CFT(j)%Freq,ctmp(j), &
       CFT(j+1)%Seq,CFT(j+1)%AA1,CFT(j+1)%Freq,ctmp(j+1), &
       CFT(j+2)%Seq,CFT(j+2)%AA1,CFT(j+2)%Freq,ctmp(j+2), &
       CFT(j+3)%Seq,CFT(j+3)%AA1,CFT(j+3)%Freq,ctmp(j+3)
    IF (MOD(i,4).eq.0)WRITE(UNIT=num,FMT="(' ')")
  END DO

! codon frequence histogram:

  DO i=1,12
    CF(i)=0
  END DO
  DO i=1,PROTlen
    tmp = CFT(CurrDNA%prot2cod(i))%Freq
    IF ((tmp.ge.0.00).and.(tmp.lt.0.05)) CF(1) = CF(1)+1
    IF ((tmp.ge.0.05).and.(tmp.lt.0.10)) CF(2) = CF(2)+1
    IF ((tmp.ge.0.10).and.(tmp.lt.0.15)) CF(3) = CF(3)+1
    IF ((tmp.ge.0.15).and.(tmp.lt.0.20)) CF(4) = CF(4)+1
    IF ((tmp.ge.0.20).and.(tmp.lt.0.25)) CF(5) = CF(5)+1
    IF ((tmp.ge.0.25).and.(tmp.lt.0.30)) CF(6) = CF(6)+1
    IF ((tmp.ge.0.30).and.(tmp.lt.0.35)) CF(7) = CF(7)+1
    IF ((tmp.ge.0.35).and.(tmp.lt.0.40)) CF(8) = CF(8)+1
    IF ((tmp.ge.0.40).and.(tmp.lt.0.45)) CF(9) = CF(9)+1
    IF ((tmp.ge.0.45).and.(tmp.lt.0.50)) CF(10) = CF(10)+1
    IF (tmp.ge.0.5) CF(11) = CF(11)+1
    CF(12) = CF(12)+1
  END DO

  WRITE(UNIT=num,FMT="(12x,' ')")
  WRITE(UNIT=num,FMT="(12x,'  Frequency Range   Number of Codons')")
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(12x,'    0% -  4%',15x,i4)") CF(1)
  WRITE(UNIT=num,FMT="(12x,'    5% -  9%',15x,i4)") CF(2)
  WRITE(UNIT=num,FMT="(12x,'   10% - 14%',15x,i4)") CF(3)
  WRITE(UNIT=num,FMT="(12x,'   15% - 19%',15x,i4)") CF(4)
  WRITE(UNIT=num,FMT="(12x,'   20% - 24%',15x,i4)") CF(5)
  WRITE(UNIT=num,FMT="(12x,'   25% - 29%',15x,i4)") CF(6)
  WRITE(UNIT=num,FMT="(12x,'   30% - 34%',15x,i4)") CF(7)
  WRITE(UNIT=num,FMT="(12x,'   35% - 39%',15x,i4)") CF(8)
  WRITE(UNIT=num,FMT="(12x,'   40% - 44%',15x,i4)") CF(9)
  WRITE(UNIT=num,FMT="(12x,'   45% - 49%',15x,i4)") CF(10)
  WRITE(UNIT=num,FMT="(12x,'      >= 50%',15x,i4)") CF(11)
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(12x,'  Total Codons Used = ',i5)") CF(12)

END SUBROUTINE Print_Codon_Histogram
SUBROUTINE Print_Codon_Log(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,n,num
  REAL :: MinFreq
  REAL :: b(64)
  CHARACTER(LEN=1) :: a(64)

  IF (TEST0) PRINT *,"Print_Codon_Log" !TEST0

! Codon Frequency Table supplied by user

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="('             CODON FREQUENCY TABLE: ',a30)") Organism
  WRITE(UNIT=num,FMT="(1x,a)") bar64
  WRITE(UNIT=num,FMT="(' ')")
  DO i=1,16
    j=((i-1)*4)+1
    WRITE(UNIT=num,FMT="(6x,a3,1x,a1,1x,f5.3,3x,a3,1x,a1,1x,f5.3,3x,a3,1x,&
      &a1,1x,f5.3,3x,a3,1x,a1,1x,f5.3)") &
                       CFT(j)%Seq,CFT(j)%AA1,CFT(j)%Freq, &
                       CFT(j+1)%Seq,CFT(j+1)%AA1,CFT(j+1)%Freq, &
                       CFT(j+2)%Seq,CFT(j+2)%AA1,CFT(j+2)%Freq, &
                       CFT(j+3)%Seq,CFT(j+3)%AA1,CFT(j+3)%Freq
    IF (MOD(i,4).eq.0)WRITE(UNIT=num,FMT="(' ')")
  END DO

! Active codons for sequence generation

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="('             ACTIVE CODONS FOR SEQUENCE GENERATION ')")
  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="('  Residue   Codons   Active Codons     Min. Freq.    Max. Freq.')")
  WRITE(UNIT=num,FMT="(1x,a)") bar64
  DO i=1,21
    MinFreq=2.0
    DO j=1,AAT(i)%NumOfCodons
      IF(CFT(AAT(i)%Codon(j))%Check) THEN
        IF(AAT(i)%Freq(j).LT.MinFreq) MinFreq=AAT(i)%Freq(j)
      END IF
    END DO
    WRITE(UNIT=num,FMT="(3x,a1,2x,a3,5x,i2,10x,i2,13x,f5.3,9x,f5.3)") &
         AAT(i)%AA1,AAT(i)%AA3,AAT(i)%NumOfCodons,AAT(i)%NumOfActiveCodons, &
         MinFreq,AAT(i)%Freq(1)
    IF (MOD(i,3).EQ.0) WRITE(UNIT=num,FMT="(' ')")
  END DO

END SUBROUTINE Print_Codon_Log
SUBROUTINE Print_Estimated_Time(SolutionNo)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER,EXTERNAL :: CurrentTimeSeconds
  INTEGER :: SolutionNo,timediff
  REAL :: current_rate,elapsed_time,time_remaining

  IF (TEST0) PRINT *,"Print_Estimated_Time" !TEST0

  timediff=CurrentTimeSeconds()-MainTimeStart
  elapsed_time=(REAL(timediff))/60
  current_rate=elapsed_time/(REAL(SolutionNo))
  time_remaining=current_rate*(REAL(TotalNumberOfSolutions-SolutionNo))
  WRITE(UNIT=console,FMT="('Elapsed time : ',f7.1,' min.')") elapsed_time
  WRITE(UNIT=console,FMT="('Finished trial',i4,1x,': rate = ',f7.1,&
    &' min/cyc : time remaining = ',f7.1,' min.')") &
    &SolutionNo,current_rate,time_remaining

  CALL FLUSH(console)

END SUBROUTINE Print_Estimated_Time
SUBROUTINE Print_FinalDNA_Log(num,SolutionNo)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j
  INTEGER :: ns
  INTEGER :: ne
  INTEGER :: nlines
  INTEGER :: SolutionNo
  INTEGER :: num
  INTEGER :: text_size
  INTEGER :: CompOligStart
  INTEGER :: half_oligos
  INTEGER :: Snext,Anext,Sprev,Aprev
  CHARACTER(LEN=9999) :: TempSeq
  CHARACTER(LEN=9999) :: mutation_array
  CHARACTER(LEN=9999) :: repeat_array
  CHARACTER(LEN=9999) :: misprime_array
  CHARACTER(LEN=9999) :: GC_array
  CHARACTER(LEN=9999) :: AT_array
  CHARACTER(LEN=9999) :: DNAlocal
  CHARACTER(LEN=9999) :: DNAcomp
  CHARACTER(LEN=9999) :: SLabel
  CHARACTER(LEN=9999) :: ALabel
  CHARACTER(LEN=9999) :: AAseq
  CHARACTER(LEN=9999) :: text,text1,text2
  LOGICAL :: repeat_present(170)     ! These are markers for each line of
  LOGICAL :: misprime_present(170)   ! output.  9999/60 = 166.65
  LOGICAL :: GC_present(170)
  LOGICAL :: AT_present(170)
  LOGICAL :: mutation_present(170)

  IF (TEST0) PRINT *,"Print_FinalDNA_Log" !TEST0

! reset global weights and parameters to defaults for final score evaluation

  MutProtPos=0
  Twt=1.0
  Cwt=1.0
  Rwt=1.0
  Mwt=1.0
  Gwt=1.0
  Awt=1.0
  Lwt=1.0
  Pwt=1.0
  Fwt=1.0

  CALL Evaluate_Scores ! one last time, to get the true value

! Store the important information

  FinalScore(SolutionNo)%FinaScore=CurrDNA%OverallScore
  FinalScore(SolutionNo)%NumOligs= CurrDNA%NumOlaps+1

! Print information for the logfile

  WRITE(UNIT=num,FMT="(' ')")
    WRITE(UNIT=num,FMT="(a80)") bar80
    WRITE(UNIT=num,FMT="(' ')")
    IF (MutantRun) THEN
      WRITE(UNIT=num,FMT="(' The DNA sequence #',i4,' (mutations in lower case):')") SolutionNo
    ELSE

! don't mess with the string, as it required for reading oldlogfiles

      WRITE(UNIT=num,FMT="(' The DNA sequence #',i4,' is:')") SolutionNo
    END IF
  WRITE(UNIT=num,FMT="(1x,a)") bar64

  nlines=DNAlen/60

! Write the DNA sequence with 60 nt per line

  TempSeq=CurrDNA%DNAseq

! Mark out mutations

  IF (MutantRun) THEN
    DO i=1,DNAlen
      IF (CurrDNA%DNAseq(i:i).ne.OLDDNAseq(i:i)) THEN
        CALL ToLowerCase(TempSeq(i:i))
      END IF
    END DO
  END IF

  DO i=1,nlines
    ne=i*60
    ns=ne-59
    WRITE(UNIT=num,FMT="(i4,1x,a60)") ns,TempSeq(ns:ne)
  END DO

! Now write the last line of the sequence (less than 60 nt)

  ns=nlines*60+1                  ! find the current nt number
  ne=DNAlen                       ! find the total length of sequence
  WRITE(UNIT=num,FMT="(i4,1x,a)") ns,TempSeq(ns:ne)
  WRITE(UNIT=num,FMT="(1x,a)") bar64

! This next section will create and write out a double stranded DNA
!   sequence with the oligonucleotides identified, written in alternating
!   upper and lower case

  CALL Create_Final_Arrays

  SLabel=''
  ALabel=''
  repeat_array=''
  mutation_array=''
  misprime_array=''
  GC_array=''
  AT_array=''

  DO i=1,90
    repeat_present(i) = .FALSE.
    misprime_present(i) = .FALSE.
    AT_present(i) = .FALSE.
    GC_present(i) = .FALSE.
    mutation_present(i) = .FALSE.
  END DO

! Test for TBIO oligos:
!
! The first half of the oligos are forward, the second half are revcompl:

  IF (TBIO) THEN
    half_oligos=(CurrDNA%NumOlaps+1)/2
    CompOligStart = half_oligos+1
  END IF

  DO i=1,DNAlen

    DNAlocal(i:i)=CurrDNA%DNAseq(i:i)
    DNAcomp(i:i)=CurrDNA%DNAseq(i:i)
    CALL ComplStr(DNAcomp(i:i))

! Make the base lowercase for every other oligo on both strands

    IF (MOD((nt2Solig(i)+1),4).eq.0) CALL ToLowerCase(DNAlocal(i:i))
    IF (MOD(nt2Aolig(i),4).eq.0) CALL ToLowerCase(DNAcomp(i:i))

! Make the gaps and overhangs blank

    IF (nt2Solig(i).eq.0) DNAlocal(i:i)=' '
    IF (nt2Aolig(i).eq.0) DNAcomp(i:i)=' '

! Fill in the amino acids

    IF (nt2aa(i).ne.0) AAseq(i:i)=AAT(nt2aa(i))%AA1
    IF (nt2aa(i).eq.0) AAseq(i:i)=' '

! Make the oligo labels

    IF (i.eq.1) THEN       ! special case: i=1
      Aprev=0
      Sprev=0
    ELSE
      Sprev=nt2Solig(i-1)
      Aprev=nt2Aolig(i-1)
    END IF

    IF (i.eq.DNAlen) THEN  ! special case: i=DNAlen
      Anext=0
      Snext=0
    ELSE
      Snext=nt2Solig(i+1)
      Anext=nt2Aolig(i+1)
    END IF

! If in TBIO mode, there are two cases

    IF (TBIO) THEN

! If the oligos are prior to the midpoint, generate complements for lower

      IF (nt2Aolig(i).le.half_oligos) THEN
        DNAcomp(i:i)=CurrDNA%DNAseq(i:i)
        IF (nt2Aolig(i).eq.0) DNAcomp(i:i)=' '
        IF (MOD(nt2Aolig(i),4).eq.0) CALL ToLowerCase(DNAcomp(i:i))
      END IF

! If the oligos are past the midpoint, generate complements for upper

      IF (nt2Solig(i).gt.half_oligos) THEN
        DNAlocal(i:i)=CurrDNA%DNAseq(i:i)
        CALL ComplStr(DNAlocal(i:i))
        IF (nt2Solig(i).eq.0) DNAlocal(i:i)=' '
        IF (MOD((nt2Solig(i)+1),4).eq.0) CALL ToLowerCase(DNAlocal(i:i))
      END IF

! If the nucleotides are prior to the mipoint,

      IF (i.lt.CurrDNA%OlapsPos(half_oligos,1)) THEN

        IF (MOD(half_oligos,2).eq.0) THEN
          IF ((nt2Solig(i).ne.0).and.(nt2Solig(i).ne.Sprev).and.&
              &(nt2Solig(i).ne.half_oligos)) THEN
            WRITE(text,FMT="(i3,' --->')") nt2Solig(i)
            text=TRIM(text)
            text=ADJUSTL(text)
            SLabel(i:(LEN_TRIM(text))+i-1)=text(1:(LEN_TRIM(text)))
          END IF
        ELSE
          IF ((nt2Solig(i).ne.0).and.(nt2Solig(i).ne.Sprev).and.&
              &(nt2Solig(i).ne.(half_oligos+1))) THEN
            WRITE(text,FMT="(i3,' --->')") nt2Solig(i)
            text=TRIM(text)
            text=ADJUSTL(text)
            SLabel(i:(LEN_TRIM(text))+i-1)=text(1:(LEN_TRIM(text)))
          END IF
        END IF

        IF ((nt2Aolig(i).ne.0).and.(nt2Aolig(i).ne.Aprev)) THEN
          WRITE(text,FMT="(i3,' --->')") nt2Aolig(i)
          text=TRIM(text)
          text=ADJUSTL(text)
          ALabel(i:(LEN_TRIM(text))+i-1)=text(1:(LEN_TRIM(text)))
        END IF

! and if they are equal to or greater than the midpoint

      ELSE

        IF ((nt2Solig(i).ne.0).and.(nt2Solig(i).ne.Snext).and.&
           &(nt2Solig(i).ne.half_oligos)) THEN
          WRITE(text,FMT="('<---',i3)") nt2Solig(i)
          text=TRIM(text)
          SLabel(i-(LEN_TRIM(text))+1:i)=text(1:(LEN_TRIM(text)))
        END IF

        IF ((nt2Aolig(i).ne.0).and.(nt2Aolig(i).ne.Anext).and.&
          &(nt2Aolig(i).ne.half_oligos)) THEN
          WRITE(text,FMT="('<---',i3)") nt2Aolig(i)
          text=TRIM(text)
          ALabel(i-(LEN_TRIM(text))+1:i)=text(1:(LEN_TRIM(text)))
        END IF

      END IF

! If not in TBIO mode, do the regular formatting

    ELSE

      IF ((nt2Solig(i).ne.0).and.(nt2Solig(i).ne.Sprev)) THEN
        WRITE(text,FMT="(i3,' --->')") nt2Solig(i)
        text=TRIM(text)
        text=ADJUSTL(text)
        SLabel(i:(LEN_TRIM(text))+i-1)=text(1:(LEN_TRIM(text)))
      END IF

      IF ((nt2Aolig(i).ne.0).and.(nt2Aolig(i).ne.Anext)) THEN
        WRITE(text,FMT="('<---',i3)") nt2Aolig(i)
        text=TRIM(text)
        ALabel(i-(LEN_TRIM(text))+1:i)=text(1:(LEN_TRIM(text)))
      END IF

    END IF
! Make the repeat, misprime, and GC/AT arrays

    j=((i-1)/60)+1

    IF (CurrDNA%RScore(i).gt.0) THEN
      repeat_present(j) = .TRUE.
      repeat_array(i:i) = '*'
    END IF

    IF (CurrDNA%MScore(i).gt.0) THEN
      misprime_present(j) = .TRUE.
      misprime_array(i:i) = '*'
    END IF

    IF (CurrDNA%GScore(i).gt.0) THEN
      GC_present(j) = .TRUE.
      GC_array(i:i) = '*'
    END IF

    IF (CurrDNA%AScore(i).gt.0) THEN
      AT_present(j) = .TRUE.
      AT_array(i:i) = '*'
    END IF

    IF (MutantRun) THEN
      IF (CurrDNA%DNAseq(i:i).ne.OLDDNAseq(i:i)) THEN
        mutation_present(j) = .TRUE.
        mutation_array(i:i) = OLDDNAseq(i:i)
      END IF
    END IF

  END DO

! Write out the DNA sequence in formatted fashion, one line at a time

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(' The oligonucleotide assembly is:')")
  WRITE(UNIT=num,FMT="(1x,a)") bar64
  WRITE(UNIT=num,FMT="(5x,'1       10        20        30        40        50        60')")
  WRITE(UNIT=num,FMT="(5x,'|        |         |         |         |         |         |')")
  WRITE(UNIT=num,FMT="(' ')")

! Write the first full lines

  DO i=1,nlines
    ne=i*60
    ns=ne-59
    WRITE(UNIT=num,FMT="(5x,a60)") SLabel(ns:ne)
    WRITE(UNIT=num,FMT="(i4,1x,a60)") ns,DNAlocal(ns:ne)
    IF (mutation_present(i)) WRITE(UNIT=num,FMT="(5x,a60,' < original')") mutation_array(ns:ne)
    IF (repeat_present(i)) WRITE(UNIT=num,FMT="(5x,a60,' < repeat')") repeat_array(ns:ne)
    IF (misprime_present(i)) WRITE(UNIT=num,FMT="(5x,a60,' < misprime')") misprime_array(ns:ne)
    IF (GC_present(i)) WRITE(UNIT=num,FMT="(5x,a60,' < GC rich')") GC_array(ns:ne)
    IF (AT_present(i)) WRITE(UNIT=num,FMT="(5x,a60,' < AT rich')") AT_array(ns:ne)
    WRITE(UNIT=num,FMT="(5x,a60)") DNAcomp(ns:ne)
    WRITE(UNIT=num,FMT="(5x,a60)") ALabel(ns:ne)
    WRITE(UNIT=num,FMT="(5x,a60)") AAseq(ns:ne)
    WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(5x,'|        |         |         |         |         |         |')")
    WRITE(UNIT=num,FMT="(' ')")
  END DO

! Now write the last line of the sequence

  ns=nlines*60+1
  ne=DNAlen
  WRITE(UNIT=num,FMT="(5x,a)") SLabel(ns:ne)
  WRITE(UNIT=num,FMT="(i4,1x,a)") ns,DNAlocal(ns:ne)
  IF (mutation_present(i)) WRITE(UNIT=num,FMT="(5x,a,' < original')") mutation_array(ns:ne)
  IF (repeat_present(i)) WRITE(UNIT=num,FMT="(5x,a,' < repeat')") repeat_array(ns:ne)
  IF (misprime_present(i)) WRITE(UNIT=num,FMT="(5x,a,' < misprime')") misprime_array(ns:ne)
  IF (GC_present(i)) WRITE(UNIT=num,FMT="(5x,a,' < GC rich')") GC_array(ns:ne)
  IF (AT_present(i)) WRITE(UNIT=num,FMT="(5x,a,' < AT rich')") AT_array(ns:ne)
  WRITE(UNIT=num,FMT="(5x,a)") DNAcomp(ns:ne)
  WRITE(UNIT=num,FMT="(5x,a)") ALabel(ns:ne)
  WRITE(UNIT=num,FMT="(5x,a)") AAseq(ns:ne)
  WRITE(UNIT=num,FMT="(1x,a)") bar64

END SUBROUTINE Print_FinalDNA_Log
SUBROUTINE Print_Final_Tally(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,num

  IF (TEST0.and.num.eq.6) PRINT *,"Print_Final_Tally" !TEST0

  WRITE(UNIT=num,FMT="(' ')")

  IF (TotalNumberOfSolutions.EQ.1) THEN
    WRITE(UNIT=num,FMT="('                          FINAL SUMMARY FOR ',i3,' SOLUTION   ')") TotalNumberOfSolutions
  ELSE
    WRITE(UNIT=num,FMT="('                          FINAL SUMMARY FOR ',i3,' SOLUTIONS  ')") TotalNumberOfSolutions
  END IF
  WRITE(UNIT=num,FMT="(a80)") bar80
  WRITE(UNIT=num,FMT="('  #    Tm   Len  |    Score   TmRange  Short    Long   #Olig  #Repeat #Misprime')")
! WRITE(UNIT=num,FMT="('123---123---123--|--123.567-----1.2-----123-----123-----123-----123')")
  WRITE(UNIT=num,FMT="(' ')")
  DO i=1,TotalNumberOfSolutions
    IF (FinalScore(i)%FinaScore.LT.900000) THEN
    WRITE(UNIT=num,FMT="(i3,2i6,2x,'|',1x,f8.3,5x,f5.1,5i8)") &
     i,&
     FinalScore(i)%MeltT,&
     FinalScore(i)%Oligo,&
     FinalScore(i)%FinaScore,&
     FinalScore(i)%TmRange, &
     FinalScore(i)%LowestOlap,&
     FinalScore(i)%LongestOligo,&
     FinalScore(i)%NumOligs,&
     FinalScore(i)%Repeats,&
     FinalScore(i)%Misprimes
    ELSE
      WRITE(UNIT=num,FMT="(i3,2i6,2x,'|',' trial abandoned ... even number of overlaps ')") &
        i,&
        FinalScore(i)%MeltT,&
        FinalScore(i)%Oligo
    END IF
  END DO
  WRITE(UNIT=num,FMT="(' ')")

END SUBROUTINE Print_Final_Tally
SUBROUTINE Print_GapFix_Histogram(num)
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,num,itmp,GF(12),GAVG,GMIN

! Gap length histogram
! Find average gap length
! Skip first gap (5' overhang) and last gap (3'overhang)

  DO i=2,CurrDNA%NumOlaps
    itmp=itmp+CurrDNA%OlapsPos(i,1)-CurrDNA%OlapsPos(i-1,2)-1
  END DO
  GAVG=INT(itmp/(CurrDNA%NumOlaps-1))

! Now generate dynamic histogram

! initialize array
  DO i=1,12
    GF(i)=0
  END DO
  GMIN = GAVG-10
  IF (GMIN.le.0) GMIN = 0

  DO i=2,CurrDNA%NumOlaps
    itmp=CurrDNA%OlapsPos(i,1)-CurrDNA%OlapsPos(i-1,2)-1
    
    IF (itmp.lt.(GMIN+2)) GF(1)=GF(1)+1
    IF ((itmp.ge.(GMIN+2)).and.(itmp.lt.(GMIN+4))) GF(2)=GF(2)+1
    IF ((itmp.ge.(GMIN+4)).and.(itmp.lt.(GMIN+6))) GF(3)=GF(3)+1
    IF ((itmp.ge.(GMIN+6)).and.(itmp.lt.(GMIN+8))) GF(4)=GF(4)+1
    IF ((itmp.ge.(GMIN+8)).and.(itmp.lt.(GMIN+10))) GF(5)=GF(5)+1
    IF ((itmp.ge.(GMIN+10)).and.(itmp.lt.(GMIN+12))) GF(6)=GF(6)+1
    IF ((itmp.ge.(GMIN+12)).and.(itmp.lt.(GMIN+14))) GF(7)=GF(7)+1
    IF ((itmp.ge.(GMIN+14)).and.(itmp.lt.(GMIN+16))) GF(8)=GF(8)+1
    IF ((itmp.ge.(GMIN+16)).and.(itmp.lt.(GMIN+18))) GF(9)=GF(9)+1
    IF ((itmp.ge.(GMIN+18)).and.(itmp.lt.(GMIN+20))) GF(10)=GF(10)+1
    IF ((itmp.ge.(GMIN+20)).and.(itmp.lt.(GMIN+22))) GF(11)=GF(11)+1
    IF (itmp.ge.(GMIN+22)) GF(12)=GF(12)+1
  END DO

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(12x,'    Gap Range     # of Gaps')")
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(19x,'<',i2,9x,i4)") (GMIN+2),GF(1)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+2),(GMIN+3),GF(2)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+4),(GMIN+5),GF(3)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+6),(GMIN+7),GF(4)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+8),(GMIN+9),GF(5)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+10),(GMIN+11),GF(6)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+12),(GMIN+13),GF(7)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+14),(GMIN+15),GF(8)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+16),(GMIN+17),GF(9)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+18),(GMIN+19),GF(10)
  WRITE(UNIT=num,FMT="(18x,i2,'-',i2,8x,i4)") (GMIN+20),(GMIN+21),GF(11)
  WRITE(UNIT=num,FMT="(19x,'>=',i2,8x,i4)") (GMIN+22),GF(12)
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')

END SUBROUTINE Print_GapFix_Histogram
SUBROUTINE Print_Help

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

! Print help instructions

  PRINT *,""
  PRINT *,"DNAWorks 3.2.4"
  PRINT *,"David Hoover"
  PRINT *,"May 04, 2017"
  PRINT *,""
  PRINT *,"DNAWorks takes as input nucleotide and/or protein sequences, codon"
  PRINT *,"information, and other variables, and attempts to optimize a synthetic"
  PRINT *,"gene.  It then outputs the gene with a variety of histograms and metrics"
  PRINT *,"for judging the probability of success for generating the gene by PCR.  It"
  PRINT *,"also outputs the oligonucleotide sequences required for PCR synthesis of"
  PRINT *,"the synthetic gene."
  PRINT *,""
  PRINT *,"COMMAND-LINE OPTIONS"
  PRINT *,"=============================================================================="
  PRINT *,""
  PRINT *,"The command line is as follows:"
  PRINT *,""
  PRINT *,"  % dnaworks [ inputfile ] [ -t0 | -t1 | -t2 | -t3 ]"
  PRINT *,""
  PRINT *,"The default inputfile is 'DNAWORKS.inp'.  All options, except for those"
  PRINT *,"on the command line, are read from the inputfile.  See below for a complete"
  PRINT *,"description of the options."
  PRINT *,""
  PRINT *,"The flags -t0, -t1, -t2, and -t3 are for testing purposes.  They report"
  PRINT *,"the internal actions within the program based on the level input."
  PRINT *,""
  PRINT *,"  -t0   Relatively simple output, only subroutine names"
  PRINT *,"  -t1   Most subroutine names reported"
  PRINT *,"  -t2   Heavy output, all subroutines, some functions"
  PRINT *,"  -t3   Way too much output, all subroutines and functions reported"
  PRINT *,""
  PRINT *,"INPUTFILE OPTIONS"
  PRINT *,"=============================================================================="
  PRINT *,""
  PRINT *,"The input is case insensitive, except for quoted strings.  Any string"
  PRINT *,"can be quoted, but it's not necessary unless the case must be preserved or"
  PRINT *,"if there are spaces or special characters (#,!).  The quotes can"
  PRINT *,"be single or double, but must begin and end around the intended"
  PRINT *,"string.  "
  PRINT *,""
  PRINT *,"Any text that follows a '#' or '!' is considered comments, and will"
  PRINT *,"be ignored."
  PRINT *,""
  PRINT *,"Options in the inputfile are of the following types:"
  PRINT *,""
  PRINT *,"  [ S ]   string"
  PRINT *,""
  PRINT *,"Strings are converted to uppercase, unless quoted (either "" or '')"
  PRINT *,""
  PRINT *,"  [ #I ]  integer number"
  PRINT *,"  [ #R ]  real number"
  PRINT *,""
  PRINT *,"Integers are, well, integers.  Real numbers can be floating point numbers"
  PRINT *,"(e.g., 12.345) or scientific notation (e.g., -12.36E+4)."
  PRINT *,""
  PRINT *,"  [ name ]  directive"
  PRINT *,""
  PRINT *,"Directives are special strings the enable or disable particular functions."
  PRINT *,"In general, only the first 4 or 5 characters are actually read, so they"
  PRINT *,"can be abbreviated."
  PRINT *,""
  PRINT *,"Directives must be placed flat against the left margin of the input file,"
  PRINT *,"otherwise they will be ignored."
  PRINT *,""
  PRINT *,"------------------------------------------------------------------------------"
  PRINT *,""
  PRINT *,"INPUT DIRECTIVES:"
  PRINT *,""
  PRINT *,"  [ tbio ]"
  PRINT *,""
  PRINT *,"The method of gene synthesis employed by DNAWorks is termed "
  PRINT *,"'thermodynamically balanced', in that all the oligonucleotides should "
  PRINT *,"assemble and anneal at the same temperature.  The amplification occurs "
  PRINT *,"everywhere at once, and ideally can generate the gene with just one round "
  PRINT *,"of PCR.  However, there are sticky cases where the gene does not amplify, "
  PRINT *,"and constructing the gene in pieces is not successful."
  PRINT *,""
  PRINT *,"A more controlled method of gene synthesis, termed 'thermodynamically "
  PRINT *,"balanced inside-out', was developed for cases where problems occurred "
  PRINT *,"during PCR synthesis (Gao, et al., 2003). In an assembly set of "
  PRINT *,"oligonucleotides, the first half of the oligos are all synthesized in the "
  PRINT *,"sense orientation, and the other half are synthesized as reverse complements "
  PRINT *,"in the anti-sense orientation of the gene. The gene assembly and amplification "
  PRINT *,"is thus done in steps of 0.4-0.6 kb from the center pair of "
  PRINT *,"oligonucleotides outward."
  PRINT *,""
  PRINT *,"Enabling tbio will enable thermodynamically balanced inside-out output."
  PRINT *,""
  PRINT *,""
  PRINT *,"  [ nogaps ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"By default, DNAWorks will try to keep all oligos the same size as the chosen"
  PRINT *,"length.  If the size is beyond the sizes required for the chosen Tm, gaps"
  PRINT *,"are introduced between overlap regions.  The directive nogaps will keep oligos"
  PRINT *,"as short as possible, with no gaps between the overlap regions."
  PRINT *,""
  PRINT *,"Restricting oligos to no gaps may slow down the optimization somewhat, and "
  PRINT *,"may result in higher scores due to a higher probability of misprimes."
  PRINT *,""
  PRINT *,"------------------------------------------------------------------------------"
  PRINT *,""
  PRINT *,"INPUT OPTIONS:"
  PRINT *,""
  PRINT *,""
  PRINT *,"  logfile [ S ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"The default output file is 'LOGFILE.txt'.  Entering a string after the"
  PRINT *,"logfile option will change the name of the logfile."
  PRINT *," "
  PRINT *," "
  PRINT *,"  title [ S ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"It's always good to give the output a title to keep it unique and to give"
  PRINT *,"you an easy way to keep track of what the output is."
  PRINT *,""
  PRINT *,""
  PRINT *,"  timelimit [ #I ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"Set a time limit for the run, in seconds.  This keeps the program from "
  PRINT *,"running forever.  A value of 0 (the default) means no limit."
  PRINT *,""
  PRINT *,""
  PRINT *,"  solutions [ #I ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"Normally DNAWorks only generates a single solution for a set of parameters."
  PRINT *,"Since the optimization involves a lot of random number calls, and that it is"
  PRINT *,"impossible to get to the 'true minimum' by Monte Carlo methods, sometimes"
  PRINT *,"generating more than one solutions is a good thing.  Look for the best"
  PRINT *,"solution in the end.  The range is 1-99."
  PRINT *,""
  PRINT *,""
  PRINT *,"  melting [ #I ] [ low #I high #I ] [ tolerance #I ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"This governs the chosen melting or annealing temperature for the oligos."
  PRINT *,"Giving a single integer (between 55 and 75) will generate a single solution."
  PRINT *,"A range of melting temperatures can be given with the low and high options,"
  PRINT *,"and a solution for each temperature will be generated.  The tolerance value"
  PRINT *,"is by default +/- 1 degree, but it can be modified.  Don't set it too high"
  PRINT *,"or the point of the program can be lost!"
  PRINT *,""
  PRINT *,""
  PRINT *,"  length [ #I ] [ low #I high #I ] [ random ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"This sets the ideal length of the oligo.  Because the oligos can have gaps,"
  PRINT *,"they can be as long as you wish, but remember that errors accumulate in"
  PRINT *,"synthetic DNA oligos very quickly beyond around 50 nts!"
  PRINT *,""
  PRINT *,"By default, an attempt is made to force all oligos to be the same size as the"
  PRINT *,"chosen length.  On occasion this can lead to a higher probability of "
  PRINT *,"misprimes.  Also, this can limit successful optimization when sequences "
  PRINT *,"are gapfixed (see below), since gap position and size will be limited.  In"
  PRINT *,"this case, enabling the length directive random causes oligos to be "
  PRINT *,"designed with random length (between 20 nt and the length chosen)."
  PRINT *,""
  PRINT *,""
  PRINT *,"  frequency [ threshold #I ] [ random ] [ strict ] [ score ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"The frequency threshold is the cutoff for which codons are used for"
  PRINT *,"reverse translation of protein sequences into DNA.  For example, a value of"
  PRINT *,"20 will allow only those codons whose frequencies equal or exceed 20%."
  PRINT *,""
  PRINT *,"By default, DNAWorks uses the highest frequency codons for the initial"
  PRINT *,"reverse translation of the protein sequences.  Having the random option"
  PRINT *,"present causes the program to choose the initial codons at random."
  PRINT *,""
  PRINT *,"By default, DNAWorks always uses the two highest frequency codons for "
  PRINT *,"optimization.  To override this default, enabling strict will "
  PRINT *,"force the program to strictly use only those codons that are within the "
  PRINT *,"chosen codon frequency threshold.  Be careful, because setting a high "
  PRINT *,"codon frequency threshold (>20%) and strict will result in many protein "
  PRINT *,"residues with a single codon available, and thus very little room for "
  PRINT *,"optimization."
  PRINT *,""
  PRINT *,"To accelerate convergence, DNAWorks does not continuously score codon "
  PRINT *,"frequency. This is allowed because only the highest frequency codons are "
  PRINT *,"usually used.  However, for the particularly picky user, enabling scored will "
  PRINT *,"force the program to continuously evaluate the codon frequency score. This "
  PRINT *,"will have the effect of increasing the overall frequency of codons (at "
  PRINT *,"the cost of other scores...). "
  PRINT *,""
  PRINT *,""
  PRINT *,"  concentration [ oligo #R ] [ sodium #R ] [ magnesium #R ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"The concentration of oligonucleotides, monovalent cations (Na+, K+), and "
  PRINT *,"magnesium in the PCR reaction can have profound effects on the annealing "
  PRINT *,"temperatures of the oligonucleotides.  The user can enter the desired "
  PRINT *,"concentrations for the PCR reaction."
  PRINT *,""
  PRINT *,"The effects of these components on the annealing temperature is based on "
  PRINT *,"the program HyTher (Nicolas Peyret, Pirro Saro and John SantaLucia, Jr.)."
  PRINT *,""
  PRINT *,"Values are in moles per liter, and can be entered in scientific notation "
  PRINT *,"for simplicity."
  PRINT *,""
  PRINT *,"Oligonucleotides must be between 100 um (1E-4 M) and 1 nm (1E-9 M), "
  PRINT *,"monovalent cations must be between 10 and 1000 mM, and magnesium must be "
  PRINT *,"between 0 and 200 mM."
  PRINT *,""
  PRINT *,""
  PRINT *,"  repeat [ #I ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"DNAWorks continuously monitors the synthetic gene for any repeats that "
  PRINT *,"occur within the gene.  A repeat can be a direct repeat, an inverted"
  PRINT *,"repeat (which can result in a hairpin), or a palindromic repeat.  If a "
  PRINT *,"repeat occurs that is above a certain length, it can lead to stable "
  PRINT *,"annealing of oligos to unexpected positions and mispriming.  Such mispriming"
  PRINT *,"can result in either no PCR product, or a long smear on a gel."
  PRINT *,""
  PRINT *,"The value for repeat governs the minimum length of nucleotides considered"
  PRINT *,"a repeat.  The default value is 8.  Increasing this number will"
  PRINT *,"decrease the number of repeats found, while decreasing it will do the"
  PRINT *,"opposite."
  PRINT *,""
  PRINT *,""
  PRINT *,"  misprime [ #I ] [ tip #I ] [ max #I ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"The major flaw to PCR-based gene synthesis is mispriming.  This occurs when"
  PRINT *,"an oligo anneals to an unexpected position on the PCR template.  To prevent"
  PRINT *,"this from happening, DNAWorks compares the ends of each oligo with the"
  PRINT *,"current synthetic sequence and analyzes its potential to anneal to that"
  PRINT *,"site.  "
  PRINT *,""
  PRINT *,"A misprime is a special variant of a repeat, in that it only occurs at the"
  PRINT *,"business end (3') of an oligo.  "
  PRINT *,""
  PRINT *,"The first number for misprime is the length of the sequence to compare.  The "
  PRINT *,"default value is 18."
  PRINT *,""
  PRINT *,"The tip number is number of nucleotides that must be exactly identical at"
  PRINT *,"the tip of the oligo.  The default is 6.  This value is based on little more"
  PRINT *,"than guessing, but increasing it will cause very few misprimes to be"
  PRINT *,"identified, and decreasing will cause too many to be identified."
  PRINT *,""
  PRINT *,"The max number is the maximum number of non-identical nucleotides in the"
  PRINT *,"misprime sequence.  The default is 8.  This number is again a guess.  It"
  PRINT *,"is generally not understood why non-identical sequences anneal to each other,"
  PRINT *,"but it is based on structural and electrostatic principles that are way too"
  PRINT *,"difficult to incorporate into this program.  Again, increasing the number"
  PRINT *,"results in too many misprimes to be identified, decreasing it causes too few."
  PRINT *,""
  PRINT *,"Needless to say, the misprime value is just plain prudence, but not "
  PRINT *,"necessarily fact."
  PRINT *,""
  PRINT *,""
  PRINT *,"  weight [ twt #R ] [ cwt #R ] [ rwt #R ] [ mwt #R ] [ gwt #R ] [ awt #R ]"
  PRINT *,"    [ lwt #R ] [ pwt #R ] [ fwt #R ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"DNAWorks optimizes a synthetic gene by evaluating the scores of a set of "
  PRINT *,"features: annealing temperature (T), codon frequency (C), repeat (R), "
  PRINT *,"misprime potential (M), GC- (G) and AT- (A) content, length (L), gapfix (F)"
  PRINT *,"and pattern constraining (P).  The default weights of each individual feature "
  PRINT *,"score are set to 1.  By increasing the weight of an individual feature, the "
  PRINT *,"final output can be nudged to favoring one feature over the others.  For "
  PRINT *,"example, in the case where the potential synthetic genes for a set of "
  PRINT *,"sequences chronically suffers from high number of repeats, increasing the "
  PRINT *,"weight of the repeat score (RWT) might decrease the final repeat score at "
  PRINT *,"the expense of the other feature scores."
  PRINT *,""
  PRINT *,"Beware, as modulating the weights is not fully tested.  Remember that this "
  PRINT *,"merely skews the results toward one feature or another, and may do more "
  PRINT *,"harm than good.  In most cases keeping the weights balanced is the best "
  PRINT *,"approach."
  PRINT *,""
  PRINT *,""
  PRINT *,"  previous [ #I ] [ S ]"
  PRINT *,""
  PRINT *,""
  PRINT *,"DNAWorks allows old sets of oligonucleotides to be read back with a new,"
  PRINT *,"mutant gene.  It then calculates scores for the mutant gene with overlap"
  PRINT *,"positions and parameters identical to the original solution.  It then"
  PRINT *,"outputs only those oligonucleotides that need to be changed.  This is very"
  PRINT *,"useful for generating mutants, since in general only one or two new oligos"
  PRINT *,"need to be synthesized."
  PRINT *,""
  PRINT *,"The integer refers to the previous solution number, and the string is the"
  PRINT *,"name of the previous logfile."
  PRINT *,""
  PRINT *,"------------------------------------------------------------------------------"
  PRINT *,"  "
  PRINT *,"INPUT SECTIONS:"
  PRINT *,""
  PRINT *,"  nucleotide [ reverse | gapfix ]"
  PRINT *,"    ..."
  PRINT *,"  //"
  PRINT *,""
  PRINT *,"Nucleotide sequences can only include A,C,G, or T in the nucleotide"
  PRINT *,"section.  They can also include degenerate sequences:"
  PRINT *,""
  PRINT *,"    B = C or G or T         rev. compl. = V"
  PRINT *,"    D = A or G or T         rev. compl. = H"
  PRINT *,"    H = A or C or T         rev. compl. = D"
  PRINT *,"    K = G or T              rev. compl. = M"
  PRINT *,"    M = A or C              rev. compl. = K"
  PRINT *,"    N = A or C or G or T    rev. compl. = N"
  PRINT *,"    R = A or G              rev. compl. = Y"
  PRINT *,"    S = C or G              rev. compl. = S"
  PRINT *,"    V = A or C or G         rev. compl. = B"
  PRINT *,"    W = A or T              rev. compl. = W"
  PRINT *,"    Y = C or T              rev. compl. = R"
  PRINT *,""
  PRINT *,"  protein [ reverse | gapfix ]"
  PRINT *,"    ..."
  PRINT *,"  //"
  PRINT *,""
  PRINT *,"Protein sequences can be input through the protein section, but can only "
  PRINT *,"include the single-letter abbreviations of the 20 standard amino acids "
  PRINT *,"(A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y).  Stop codons are designated by X."
  PRINT *,""
  PRINT *,"The reverse directive causes the nucleotide sequence (either original or"
  PRINT *,"translated from the protein sequence) to be reversed on incorporation in"
  PRINT *,"the synthetic gene."
  PRINT *,""
  PRINT *,"The gapfix directive is used when the sequence should not fall within "
  PRINT *,"overlap regions, but rather only in the gaps or overhangs that are single"
  PRINT *,"stranded in the annealed assembly prior to PCR.  This is advantageous for"
  PRINT *,"subsequent mutations by oligonucleotide replacement. For example, if a "
  PRINT *,"synthetic gene will be exhaustively mutated at a single codon, having the"
  PRINT *,"codon entirely within a gap region will allow its mutation by replacing a"
  PRINT *,"single oligonucleotide, rather than two or three."
  PRINT *,""
  PRINT *,"The gapfix directive will enable Fixed Gap Scoring.  Any nt that are"
  PRINT *,"designated as gapfixed but fall within overlap regions will increase the"
  PRINT *,"global score.  DNAWorks will then try to minimize the score by moving the"
  PRINT *,"gap regions toward the gapfixed nucleotides.  Because gap regions are "
  PRINT *,"generally short (less than 10 nt), the sequence should be very short.  "
  PRINT *,"Otherwise the global score will remain quite high, and other features (Tm,"
  PRINT *,"repeats, misprimes) will not receive as much attention."
  PRINT *,""
  PRINT *,"Gapfixing is much more effective when oligo lengths are allowed be "
  PRINT *,"randomized, rather than fixed to the length chosen by default.  See"
  PRINT *,"length option, above, for more details."
  PRINT *,""
  PRINT *,""
  PRINT *,"  codon [ ecoli2 | E. coli | C. elegans | D. melanogaster | H. sapiens | "
  PRINT *,"    M. musculus | R. novegicus | S. cerevesiae | X. laevis | P. pastoris ]"
  PRINT *,"    [ ..."
  PRINT *,"  // ]"
  PRINT *,""
  PRINT *,"Codon frequencies can be entered manually in the codon section using "
  PRINT *,"GCG-format codon frequencies.  If a directive corresponding to a given "
  PRINT *,"organism is present, the codon frequency for that organism will be used."
  PRINT *,""
  PRINT *,"  pattern"
  PRINT *,"    ..."
  PRINT *,"  //"
  PRINT *,""
  PRINT *,"Nucleotide patterns can be screened if entered in the pattern section."
  PRINT *,"Pattern sequences can be normal or degenerate nucleotide sequences."
  PRINT *,""

  STOP

END SUBROUTINE Print_Help
SUBROUTINE Print_Histogram(num,SolutionNo)
!
! This subroutine prints out four histograms, a codon use histogram, a melting
! temperature histogram, a length histogram, and a repeats histogram.
! These are useful for determining the real quality of a solution.
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: num,SolutionNo

  IF (TEST0) PRINT *,"Print_Histogram" !TEST0
  IF (PROTlen.gt.0) CALL Print_Codon_Histogram(num)
  CALL Print_Melt_Histogram(num,SolutionNo)
  CALL Print_Overlap_Histogram(num,SolutionNo)
  CALL Print_Length_Histogram(num,SolutionNo)
  IF (GapFix) CALL Print_GapFix_Histogram(num)

  CALL Export_Misprime_Pairs(num)
  CALL Export_Repeat_Pairs(num)

  FinalScore(SolutionNo)%Repeats = CurrDNA%RN
  FinalScore(SolutionNo)%Misprimes = CurrDNA%MSN

  WRITE(UNIT=num,FMT="(' ')")

END SUBROUTINE Print_Histogram
SUBROUTINE Print_Length_Histogram(num,SolutionNo)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,OligoSize(999),num,itmp,SolutionNo,NT(12),NTLONG,NAVG,NMIN

  IF (TEST0) PRINT *,"Print_Length_Histogram" !TEST0

! Dynamic Length histogram

  NTLONG=0

! The first oligo

  OligoSize(1)=CurrDNA%OlapsPos(1,2)
  itmp=OligoSize(1)

  j=2
  DO i=3,CurrDNA%NumOlaps,2      ! First find oligos in the sense strand
    OligoSize(j)=CurrDNA%OlapsPos(i,2)-CurrDNA%OlapsPos(i-1,1)+1
    itmp=itmp+OligoSize(j)
    j=j+1
  END DO

  DO i=CurrDNA%NumOlaps-1,2,-2   ! Then find oligos in the antisense strand
    OligoSize(j)=CurrDNA%OlapsPos(i,2)-CurrDNA%OlapsPos(i-1,1)+1
    itmp=itmp+OligoSize(j)
    j=j+1
  END DO

! The last oligo

  OligoSize(j)=DNAlen-CurrDNA%OlapsPos(CurrDNA%NumOlaps,1)+1
  itmp=itmp+OligoSize(j)

  NAVG=INT(itmp/(CurrDNA%NumOlaps+1))

!     initialize array
  DO i=1,12
    NT(i)=0
  END DO
  NMIN = NAVG-10
  IF (NMIN.le.0) NMIN = 0

  DO i=1,CurrDNA%NumOlaps+1
    itmp=OligoSize(i)
    IF (itmp.gt.NTLONG) NTLONG=itmp
    IF (itmp.lt.(NMIN+2)) NT(1)=NT(1)+1
    IF ((itmp.ge.(NMIN+2)).and.(itmp.lt.(NMIN+4))) NT(2)=NT(2)+1
    IF ((itmp.ge.(NMIN+4)).and.(itmp.lt.(NMIN+6))) NT(3)=NT(3)+1
    IF ((itmp.ge.(NMIN+6)).and.(itmp.lt.(NMIN+8))) NT(4)=NT(4)+1
    IF ((itmp.ge.(NMIN+8)).and.(itmp.lt.(NMIN+10))) NT(5)=NT(5)+1
    IF ((itmp.ge.(NMIN+10)).and.(itmp.lt.(NMIN+12))) NT(6)=NT(6)+1
    IF ((itmp.ge.(NMIN+12)).and.(itmp.lt.(NMIN+14))) NT(7)=NT(7)+1
    IF ((itmp.ge.(NMIN+14)).and.(itmp.lt.(NMIN+16))) NT(8)=NT(8)+1
    IF ((itmp.ge.(NMIN+16)).and.(itmp.lt.(NMIN+18))) NT(9)=NT(9)+1
    IF ((itmp.ge.(NMIN+18)).and.(itmp.lt.(NMIN+20))) NT(10)=NT(10)+1
    IF ((itmp.ge.(NMIN+20)).and.(itmp.lt.(NMIN+22))) NT(11)=NT(11)+1
    IF (itmp.ge.(NMIN+22)) NT(12)=NT(12)+1
  END DO

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(12x,'  Length Range   # of Oligos')")
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(19x,'<',i3,9x,i4)") (NMIN+2),NT(1)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+2),(NMIN+3),NT(2)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+4),(NMIN+5),NT(3)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+6),(NMIN+7),NT(4)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+8),(NMIN+9),NT(5)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+10),(NMIN+11),NT(6)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+12),(NMIN+13),NT(7)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+14),(NMIN+15),NT(8)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+16),(NMIN+17),NT(9)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+18),(NMIN+19),NT(10)
  WRITE(UNIT=num,FMT="(18x,i3,'-',i3,8x,i4)") (NMIN+20),(NMIN+21),NT(11)
  WRITE(UNIT=num,FMT="(19x,'>=',i3,8x,i4)") (NMIN+22),NT(12)
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(12x,'  Longest = ',i3)") NTLONG
  FinalScore(SolutionNo)%LongestOligo = NTLONG

END SUBROUTINE Print_Length_Histogram
SUBROUTINE Print_Melt_Histogram(num,SolutionNo)
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,num,SolutionNo,MT(10),MMIN
  REAL :: tmp,MAVG,MLOW,MHIGH

  IF (TEST0) PRINT *,"Print_Melt_Histogram" !TEST0

! Dynamic Melting temperature histogram

!   Find average overlap length

  tmp=0.0
  DO i=1,CurrDNA%NumOlaps
    tmp=tmp+CurrDNA%MeltT(i)
  END DO
  MAVG=tmp/CurrDNA%NumOlaps
!     initialize array

  DO i=1,10
    MT(i)=0
  END DO
  MMIN = INT(MAVG)-3
  IF (MMIN.le.0) MMIN = 0
  MHIGH=0.0
  MLOW=999.0

!   generate array

  DO i=1,CurrDNA%NumOlaps
    tmp=CurrDNA%MeltT(i)
    IF (tmp.gt.MHIGH) MHIGH=tmp
    IF (tmp.lt.MLOW) MLOW=tmp
    IF (tmp.lt.(MMIN+1)) MT(1)=MT(1)+1
    IF ((tmp.ge.(MMIN+1)).and.(tmp.lt.(MMIN+2))) MT(2)=MT(2)+1
    IF ((tmp.ge.(MMIN+2)).and.(tmp.lt.(MMIN+3))) MT(3)=MT(3)+1
    IF ((tmp.ge.(MMIN+3)).and.(tmp.lt.(MMIN+4))) MT(4)=MT(4)+1
    IF ((tmp.ge.(MMIN+4)).and.(tmp.lt.(MMIN+5))) MT(5)=MT(5)+1
    IF ((tmp.ge.(MMIN+5)).and.(tmp.lt.(MMIN+6))) MT(6)=MT(6)+1
    IF ((tmp.ge.(MMIN+6)).and.(tmp.lt.(MMIN+7))) MT(7)=MT(7)+1
    IF ((tmp.ge.(MMIN+7)).and.(tmp.lt.(MMIN+8))) MT(8)=MT(8)+1
    IF ((tmp.ge.(MMIN+8)).and.(tmp.lt.(MMIN+9))) MT(9)=MT(9)+1
    IF (tmp.ge.(MMIN+9)) MT(10)=MT(10)+1
  END DO

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(12x,'      Tm Range    # of Overlaps ')")
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(19x,'<',i2,9x,i4)") (MMIN+1),MT(1)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+1),MT(2)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+2),MT(3)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+3),MT(4)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+4),MT(5)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+5),MT(6)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+6),MT(7)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+7),MT(8)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (MMIN+8),MT(9)
  WRITE(UNIT=num,FMT="(18x,'>=',i2,9x,i4)") (MMIN+9),MT(10)
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(12x'  Tm Range = ',f4.1)") (MHIGH-MLOW)
  FinalScore(SolutionNo)%TmRange = (MHIGH-MLOW)

END SUBROUTINE Print_Melt_Histogram
SUBROUTINE Print_Oligo_Log(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: last,i,ns,ne,S,num
  INTEGER :: OligoSize(999)
  CHARACTER(LEN=9999) :: OldOligo(999)
  CHARACTER(LEN=9999) :: NewOligo(999)

  last=CurrDNA%NumOlaps+1

  IF (TEST0) PRINT *,"Print_Oligo_Log" !TEST0

! Find OligoSizes:

  OligoSize(1)=CurrDNA%OlapsPos(1,2)
  DO i=2,CurrDNA%NumOlaps
    OligoSize(i)=CurrDNA%OlapsPos(i,2)-CurrDNA%OlapsPos(i-1,1)+1
  END DO
  OligoSize(last)=DNAlen-CurrDNA%OlapsPos(CurrDNA%NumOlaps,1)+1

! Now create the oligos, starting with the first

  S=CurrDNA%OlapsPos(1,2)
  NewOligo(1)(1:S)=CurrDNA%DNAseq(1:S)
  OldOligo(1)(1:S)=OLDDNAseq(1:S)

! the middle oligos,

  DO i=2,CurrDNA%NumOlaps
    ns=CurrDNA%OlapsPos(i-1,1)
    ne=CurrDNA%OlapsPos(i,2)
    S=ne-ns+1

! However, if in TBIO mode, then

    IF (TBIO) THEN

!   write all oligos before the midpoint in sense

      IF (i.lt.((CurrDNA%NumOlaps+1)/2)+1) THEN
        NewOligo(i)(1:S) = CurrDNA%DNAseq(ns:ne)
        OldOligo(i)(1:S) = OLDDNAseq(ns:ne)

!   write later oligos in reverse complements

      ELSE
        NewOligo(i)(1:S) = CurrDNA%DNAseq(ns:ne)
        CALL RevComplStr(NewOligo(i)(1:S))
        OldOligo(i)(1:S) = OLDDNAseq(ns:ne)
        CALL RevComplStr(OldOligo(i)(1:S))
      END IF

    ELSE

      IF (MOD(i,2).eq.1) THEN
        NewOligo(i)(1:S) = CurrDNA%DNAseq(ns:ne)
        OldOligo(i)(1:S) = OLDDNAseq(ns:ne)
      ELSE
        NewOligo(i)(1:S) = CurrDNA%DNAseq(ns:ne)
        CALL RevComplStr(NewOligo(i)(1:S))
        OldOligo(i)(1:S) = OLDDNAseq(ns:ne)
        CALL RevComplStr(OldOligo(i)(1:S))
      END IF

    END IF
  END DO

! and the last oligo.

  ns=CurrDNA%OlapsPos(CurrDNA%NumOlaps,1)
  ne=DNAlen
  S=ne-ns+1
  NewOligo(last)(1:S) = CurrDNA%DNAseq(ns:ne)
  CALL RevComplStr(NewOligo(last)(1:S))
  OldOligo(last)(1:S) = OLDDNAseq(ns:ne)
  CALL RevComplStr(OldOligo(last)(1:S))

! If this is a MutantRun, only print the different oligos

  S=0
  DO i=1,last
    IF (MutantRun) THEN
      IF (NewOligo(i)(1:OligoSize(i)).ne.OldOligo(i)(1:OligoSize(i))) THEN
        S=S+1
      END IF
    END IF
  END DO

  IF (S.eq.1) THEN
    WRITE(UNIT=num,FMT="(5x,'  1 replacement oligonucleotide needs to be synthesized')")
  ELSE IF (S.gt.1) THEN
    WRITE(UNIT=num,FMT="(5x,i3,' replacement oligonucleotides need to be synthesized')") S
  ELSE
    IF (.not.MutantRun) THEN
      WRITE(UNIT=num,FMT="(7x,i3,' oligonucleotides need to be synthesized')") CurrDNA%NumOlaps+1
    ELSE
      WRITE(UNIT=num,FMT="(7x,'No new oligonucleotides need to be synthesized')")
    END IF
  END IF

  WRITE(UNIT=num,FMT="(1x,a64)") bar64

  DO i=1,last
    IF (MutantRun) THEN
      IF (NewOligo(i)(1:OligoSize(i)).ne.OldOligo(i)(1:OligoSize(i))) THEN
        WRITE(UNIT=num,FMT="(i3,1x,a,1x,i3)") i,NewOligo(i)(1:OligoSize(i)),OligoSize(i)
      END IF
    ELSE
      WRITE(UNIT=num,FMT="(i3,1x,a,1x,i3)") i,NewOligo(i)(1:OligoSize(i)),OligoSize(i)
    END IF
  END DO

END SUBROUTINE Print_Oligo_Log
SUBROUTINE Print_Output_End(num)
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: num,i,version
  CHARACTER(LEN=80),EXTERNAL :: CenterStr
  CHARACTER(LEN=80) :: a

  IF (TEST0.and.num.eq.6) PRINT *,"Print_Output_End" !TEST0

  IF (JACEK) WRITE(num,FMT='(a)') bar80
  IF (JACEK) a='Macromolecular Crystallography Laboratory,'
  IF (JACEK) WRITE(num,FMT='(a)') CenterStr(a)
  IF (JACEK) a='Molecular Assembly and Cell Signaling Section'
  IF (JACEK) WRITE(num,FMT='(a)') CenterStr(a)
  IF (JACEK) a='in collaboration with'
  IF (JACEK) WRITE(num,FMT='(a)') CenterStr(a)
  WRITE(num,FMT='(a)') bar80
  a='HPC @ NIH -- Center for Information Technology'
  WRITE(num,FMT='(a)') CenterStr(a)
  a='https://hpc.nih.gov'
  WRITE(num,FMT='(a)') CenterStr(a)
  a='National Institutes of Health, Department of Health and Human Services'
  WRITE(num,FMT='(a)') CenterStr(a)
  a=''
  WRITE(num,FMT='(a)') CenterStr(a)
  a='DNAWorks Web Site: https://hpcwebapps.cit.nih.gov/dnaworks'
  WRITE(num,FMT='(a)') CenterStr(a)
  WRITE(num,FMT='(a)') bar80

  WRITE(num,FMT='(a)') ''

END SUBROUTINE Print_Output_End
SUBROUTINE Print_Output_Start(num)
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: num
  CHARACTER(LEN=8) :: CurrentTime,x
  CHARACTER(LEN=10) :: CurrentDate,y
  CHARACTER(LEN=80),EXTERNAL :: CenterStr
  CHARACTER(LEN=80) :: a

  IF (TEST0.and.num.eq.6) PRINT *,"Print_Output_Start" !TEST0

  x=CurrentTime()
  y=CurrentDate()

  IF (JACEK) WRITE(num,FMT='(a)') bar80
  IF (JACEK) a='Macromolecular Crystallography Laboratory,'
  IF (JACEK) WRITE(num,FMT='(a)') CenterStr(a)
  IF (JACEK) a='Molecular Assembly and Cell Signaling Section'
  IF (JACEK) WRITE(num,FMT='(a)') CenterStr(a)
  IF (JACEK) WRITE(num,FMT='(a)') bar80
  IF (JACEK) a='in collaboration with'
  IF (JACEK) WRITE(num,FMT='(a)') CenterStr(a)
  WRITE(num,FMT='(a)') bar80
  a='HPC @ NIH -- Center for Information Technology'
  WRITE(num,FMT='(a)') CenterStr(a)
  a='https://hpc.nih.gov'
  WRITE(num,FMT='(a)') CenterStr(a)
  a='National Institutes of Health, Department of Health and Human Services'
  WRITE(num,FMT='(a)') CenterStr(a)
  a=''
  WRITE(num,FMT='(a)') CenterStr(a)
  a='DNAWorks Web Site: https://hpcwebapps.cit.nih.gov/dnaworks'
  WRITE(num,FMT='(a)') CenterStr(a)
  WRITE(num,FMT='(a)') bar80
  a=''
  WRITE(num,FMT='(a)') CenterStr(a)
  a='Send all correspondence to webtools@hpc.nih.gov'
  WRITE(num,FMT='(a)') CenterStr(a)
  WRITE(num,FMT='(a)') bar80

  WRITE(num,FMT='(a)') ''
  WRITE(num,FMT="(' Job started on ',a10,' at ',a8)") y,x

! jobname

  IF (LEN_TRIM(jobname).GT.0) THEN
    WRITE(num,FMT='(a)') ''
    WRITE(num,FMT="(' Job name: ',a69)") jobname
  END IF

! email

  IF (LEN_TRIM(email).GT.5) THEN
    WRITE(num,FMT='(a)') ''
    WRITE(num,FMT="(' Output will be sent to ',a30)") email
  END IF

END SUBROUTINE Print_Output_Start
SUBROUTINE Print_Overlap_Histogram(num,SolutionNo)
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,num,itmp,SolutionNo,OL(12),OAVG,OMIN,OLOW 

  IF (TEST0) PRINT *,"Print_Overlap_Histogram" !TEST0

! Dynamic overlap length histogram

!   Find average overlap length

  itmp=0
  DO i=1,CurrDNA%NumOlaps
    itmp=itmp+CurrDNA%OlapsPos(i,2)-CurrDNA%OlapsPos(i,1)+1
  END DO
  OAVG=INT(itmp/(CurrDNA%NumOlaps))

!     initialize array

  OLOW = 99
  DO i=1,12
    OL(i)=0
  END DO
  OMIN = OAVG-4
  IF (OMIN.le.0) OMIN = 0

!   generate OL array

  DO i=1,CurrDNA%NumOlaps
    itmp=CurrDNA%OlapsPos(i,2)-CurrDNA%OlapsPos(i,1)+1
    IF (itmp.lt.OLOW) OLOW = itmp

    IF (itmp.lt.(OMIN+1)) OL(1)=OL(1)+1
    IF (itmp.eq.(OMIN+1)) OL(2)=OL(2)+1
    IF (itmp.eq.(OMIN+2)) OL(3)=OL(3)+1
    IF (itmp.eq.(OMIN+3)) OL(4)=OL(4)+1
    IF (itmp.eq.(OMIN+4)) OL(5)=OL(5)+1
    IF (itmp.eq.(OMIN+5)) OL(6)=OL(6)+1
    IF (itmp.eq.(OMIN+6)) OL(7)=OL(7)+1
    IF (itmp.eq.(OMIN+7)) OL(8)=OL(8)+1
    IF (itmp.eq.(OMIN+8)) OL(9)=OL(9)+1
    IF (itmp.eq.(OMIN+9)) OL(10)=OL(10)+1
    IF (itmp.eq.(OMIN+10)) OL(11)=OL(11)+1
    IF (itmp.eq.(OMIN+11)) OL(12)=OL(12)+1
  END DO
   
  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(12x,'  Ovrlap Len Range  # of Oligos')")
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(19x,'<',i2,9x,i4)") (OMIN+1),OL(1)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+1),OL(2)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+2),OL(3)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+3),OL(4)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+4),OL(5)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+5),OL(6)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+6),OL(7)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+7),OL(8)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+8),OL(9)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+9),OL(10)
  WRITE(UNIT=num,FMT="(20x,i2,9x,i4)") (OMIN+10),OL(11)
  WRITE(UNIT=num,FMT="(18x,'>=',i2,9x,i4)") (OMIN+11),OL(12)
  WRITE(UNIT=num,FMT='(12x," -------------------------------------")')
  WRITE(UNIT=num,FMT="(12x,'  Lowest Overlap = ',i2)") OLOW
  FinalScore(SolutionNo)%LowestOlap = OLOW

END SUBROUTINE Print_Overlap_Histogram
SUBROUTINE Print_Param_Log(num,SolutionNo)
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: num,SolutionNo,i
  CHARACTER(LEN=30) :: temp

  IF (TEST0) PRINT *,"Print_Param_Log" !TEST0

! For the logfile

  WRITE(num,FMT="(' ')")
  WRITE(num,FMT="(1x,a64)") bar64

  WRITE(temp,FMT="(i4)") SolutionNo
  temp=TRIM(temp)
  temp=ADJUSTL(temp)

! don't mess with these strings, as they are required for reading oldlogfiles

  WRITE(num,FMT='(20x,"PARAMETERS FOR TRIAL ",a)') temp

  WRITE(num,FMT="(1x,a64)") bar64

  WRITE(temp,FMT="(i4,' nt')") DNAlen
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Total Size Of Gene ......... ",a)') temp

  WRITE(temp,FMT="(i4)") PROTlen
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Protein Residues ........... ",a)') temp

  i=mutPROTnum
  WRITE(temp,FMT="(i4)") i
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Mutatable Residues ......... ",a)') temp

  i=DNAlen-(mutPROTnum*3)
  WRITE(temp,FMT="(i4,' nt')") i
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Fixed Nucleotides .......... ",a)') temp

  WRITE(temp,FMT="(i4,' nt')") NumDegPos
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Degenerate Nucleotides ..... ",a)') temp

  WRITE(temp,FMT="(i4,' nt')") OligoLen
  IF (OligoLenRandom) temp=temp(1:(LEN_TRIM(temp)))//' RANDOM'
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Oligo Size ................. ",a)') temp

  WRITE(temp,FMT="(i4,' +/-',i2,'*C')") MeltTemp,MeltTol
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Annealing Temp ............. ",a)') temp

  WRITE(temp,FMT="(es8.2e1,' M')") OligoConc
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Oligo Concentration ........ ",a)') temp

  WRITE(temp,FMT="(es8.2e1,' M')") SodiumConc
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Sodium Concentration ....... ",a)') temp

  WRITE(temp,FMT="(es8.2e1,' M')") MgConc
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Mg2+ Concentration ......... ",a)') temp

  WRITE(temp,FMT="(i4,'%')") SeqOptimToler
  IF (CodonStrict) temp=temp(1:(LEN_TRIM(temp)))//' STRICT'
  IF (CodonRandom) temp=temp(1:(LEN_TRIM(temp)))//' RANDOM'
  IF (ScoreCodons) temp=temp(1:(LEN_TRIM(temp)))//' SCORED'
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Codon Frequency Threshold .. ",a)') temp

  WRITE(temp,FMT="(i2,' nt')") RepLen
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Repeat Threshold ........... ",a)') temp

  WRITE(temp,FMT="(i2,'/',i2,' (',i1,' exact) nt')") MaxNonId,MPLn,MPTip
  temp=TRIM(temp)
  temp=ADJUSTL(temp)
  WRITE(num,FMT='(15x,"Mispriming Threshold ....... ",a)') temp

  IF (Awt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Awt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"AT Score Weight ............ ",a)') temp
  END IF

  IF (Cwt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Cwt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"Codon Score Weight ......... ",a)') temp
  END IF

  IF (Fwt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Fwt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"Fixed Gap Score Weight ..... ",a)') temp
  END IF

  IF (Gwt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Gwt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"GC Score Weight ............ ",a)') temp
  END IF

  IF (Lwt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Lwt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"Length Score Weight ........ ",a)') temp
  END IF

  IF (Mwt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Mwt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"Misprime Score Weight ...... ",a)') temp
  END IF

  IF (Pwt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Pwt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"Pattern Score Weight ....... ",a)') temp
  END IF

  IF (Rwt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Rwt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"Repeat Score Weight ........ ",a)') temp
  END IF

  IF (Twt.ne.1.0) THEN
    WRITE(temp,FMT="(es8.2e1)") Twt
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    WRITE(num,FMT='(15x,"Tm Score Weight ............ ",a)') temp
  END IF

  IF (TBIO) &
    WRITE(num,FMT='(15x,"Thermodynamically Balanced Inside-Out mode output")')

  IF (NOGAPS) &
    WRITE(num,FMT='(15x,"No gaps are allowed in assembled gene")')

END SUBROUTINE Print_Param_Log
SUBROUTINE Print_Pattern_Log(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,num
  CHARACTER(LEN=30) :: text

  IF (TEST0) PRINT *,"Print_Pattern_Log" !TEST0

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT='(16x,"SEQUENCE PATTERNS TO BE SCREENED")')
  WRITE(UNIT=num,FMT="(1x,a64)") bar64

  IF (PTNnum.gt.0) THEN
    DO i=1,PTNnum
      WRITE(UNIT=num,FMT="(2x,a12,a66)") PTN(i)%Name(1:12),PTN(i)%Seq(1:66)
    END DO
  END IF
  WRITE(UNIT=num,FMT="(1x,a64)") bar64

END SUBROUTINE Print_Pattern_Log
SUBROUTINE Print_Pattern_Screen(num)
!
! This subroutine will perform a final restriction digest of the DNA using the input sites
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,n,num,SolutionNo,ln,NumHits
  CHARACTER(LEN=9999) :: text
  CHARACTER(LEN=9999) :: ftext,rtext
  LOGICAL,EXTERNAL :: DegCmpr
  INTEGER :: curr,start

  TYPE junk
    CHARACTER(LEN=80) :: S
    INTEGER :: L
    CHARACTER(LEN=80) :: N
    INTEGER :: P
    CHARACTER(LEN=40) :: T
    LOGICAL :: D
  END TYPE

  TYPE(junk) :: x(999),c

  IF (TEST0) PRINT *,"Print_Pattern_Screen" !TEST0

! Write the header

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(10x,'Sequence Patterns Screened (As Supplied By User)')")
  WRITE(UNIT=num,FMT="(1x,a)") bar64
!
! Find the sites
!
  NumHits = 0

  main: DO i=1,PTNnum
    IF (PTN(i)%Degen) THEN
      ln=PTN(i)%Len
      ftext=PTN(i)%Seq(1:PTN(i)%Len)
      rtext=PTN(i)%SeqRC(1:PTN(i)%Len)
      IF (ln.EQ.0) EXIT main
      DO j=1,(DNAlen-ln)
        IF (DegCmpr(ftext,CurrDNA%DNAseq(j:j+ln-1))) THEN
          NumHits=NumHits+1
          x(NumHits)%S=PTN(i)%Seq
          x(NumHits)%L=PTN(i)%Len
          x(NumHits)%N=PTN(i)%Name
          x(NumHits)%P=j
          x(NumHits)%T="forward"
          x(NumHits)%D=.TRUE.
        END IF
        IF (.not.PTN(i)%SelfCompl) THEN
          IF (DegCmpr(rtext,CurrDNA%DNAseq(j:j+ln-1))) THEN
            NumHits=NumHits+1
            x(NumHits)%S=PTN(i)%Seq
            x(NumHits)%L=PTN(i)%Len
            x(NumHits)%N=PTN(i)%Name
            x(NumHits)%P=j
            x(NumHits)%T="reverse"
            x(NumHits)%D=.TRUE.
          END IF
        END IF
      END DO
    ELSE  ! no degenerates are present
      curr=0
      start=1
! forward direction
      forward: DO n=1,DNAlen
        j=INDEX(CurrDNA%DNAseq(start:DNAlen),PTN(i)%Seq(1:PTN(i)%Len))
        curr=curr+j
        IF (j.eq.0) THEN
          EXIT forward
        ELSE
          NumHits=NumHits+1
          x(NumHits)%S=PTN(i)%Seq
          x(NumHits)%L=PTN(i)%Len
          x(NumHits)%N=PTN(i)%Name
          x(NumHits)%P=curr
          x(NumHits)%T="forward"
          x(NumHits)%D=.FALSE.
          start=curr+1
        END IF
      END DO forward
! reverse direction if needed
      IF (.not.PTN(i)%SelfCompl) THEN
        curr=0
        start=1
        reverse: DO n=1,DNAlen
          j=INDEX(CurrDNA%DNAseq(start:DNAlen),PTN(i)%SeqRC(1:PTN(i)%Len))
          curr=curr+j
          IF (j.eq.0) THEN
            EXIT reverse
          ELSE
            NumHits=NumHits+1
            x(NumHits)%S=PTN(i)%Seq
            x(NumHits)%L=PTN(i)%Len
            x(NumHits)%N=PTN(i)%Name
            x(NumHits)%P=curr+PTN(i)%Len-1
            x(NumHits)%T="reverse"
            x(NumHits)%D=.FALSE.
            start=curr+1
          END IF
        END DO reverse
      END IF
    END IF
  END DO main

! Sort by position

  IF (NumHits.gt.1) THEN
    DO i=(NumHits-1),1,-1                     ! type(junk) swap
      DO j=1,i
        IF (x(j)%P.gt.x((j+1))%P) THEN
          c=x(j)
          x(j)=x(j+1)
          x(j+1)=c
        END IF
      END DO
    END DO
  END IF

! Write out the sites in order

   IF (NumHits.gt.0) THEN
     WRITE(UNIT=num,FMT="(2x,'Name        Seq      Pos    Notes')")
     WRITE(UNIT=num,FMT="(1x,a)") bar64
     DO i=1,NumHits
       IF (x(i)%D) THEN
         IF (x(i)%T(8:27).eq.'') THEN
           x(i)%T(8:19) = ", degenerate"
         ELSE
           x(i)%T(28:39) = ", degenerate"
         END IF
       END IF
       WRITE(UNIT=num,FMT="(2x,a12,a9,i5,2x,a40)") x(i)%N,x(i)%S,x(i)%P,x(i)%T
     END DO
   ELSE
     WRITE(UNIT=num,FMT="(23x,'None found')")
   END IF

  WRITE(UNIT=num,FMT="(1x,a)") bar64
  WRITE(UNIT=num,FMT="(' ')")

END SUBROUTINE Print_Pattern_Screen
SUBROUTINE Print_Scores_Log(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: num

  IF (TEST0) PRINT *,"Print_Scores_Log" !TEST0

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT="(1x,&
    &'The total codon usage score ........... ',f9.3)")CurrDNA%TotalCScore
  WRITE(UNIT=num,FMT="(1x,&
    &'The total length score ................ ',f9.3)")CurrDNA%TotalLScore
  WRITE(UNIT=num,FMT="(1x,&
    &'The total melting temperature score ... ',f9.3)")CurrDNA%TotalTScore
  WRITE(UNIT=num,FMT="(1x,&
    &'The total repeat score ................ ',f9.3)")CurrDNA%TotalRScore
  WRITE(UNIT=num,FMT="(1x,&
    &'The total pattern score ............... ',f9.3)")CurrDNA%TotalPScore
  WRITE(UNIT=num,FMT="(1x,&
    &'The total mispriming score ............ ',f9.3)")CurrDNA%TotalMScore
  WRITE(UNIT=num,FMT="(1x,&
    &'The total AT content score ............ ',f9.3)")CurrDNA%TotalAScore
  WRITE(UNIT=num,FMT="(1x,&
    &'The total GC content score ............ ',f9.3)")CurrDNA%TotalGScore
  IF (GapFix) THEN
    WRITE(UNIT=num,FMT="(1x,&
      &'The total fixed gap score ............. ',f9.3)")CurrDNA%TotalFScore
  ELSE
    WRITE(UNIT=num,FMT="(1x,&
      &'The total fixed gap score ............. ','     N/A')")
  END IF
  WRITE(UNIT=num,FMT="(1x,&
    &'               The OVERALL score ...... ',f9.3)")CurrDNA%OverallScore
  WRITE(UNIT=num,FMT="(' ')")

END SUBROUTINE Print_Scores_Log
SUBROUTINE Print_Seq_Log(num)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,ns,ne,nlines,num,j,cnt,k,start,finish
  CHARACTER(LEN=30) :: temp

  IF (TEST0) PRINT *,"Print_Seq_Log" !TEST0

! For the logfile
!
  DO j=1,NumberOfSeq 

    temp = ''
    IF (SeqReverse(j)) temp=temp(1:(LEN_TRIM(temp)))//' REVERSED'
    IF (SeqGapFix(j)) temp=temp(1:(LEN_TRIM(temp)))//' GAPFIXED'
    temp=TRIM(temp)
    temp=ADJUSTL(temp)
    cnt=0
    start=0
    DO k=1,INITlen
      IF (INIT2seq(k).eq.j) THEN 
        cnt=cnt+1
        IF (start.eq.0) start=k
      END IF
    END DO
    finish=start+cnt-1
    WRITE(UNIT=num,FMT="(' ')")

    IF (SeqIsProt(j)) THEN
      WRITE(UNIT=num,FMT='(6x,"SEQUENCE ",i2,": PROTEIN LENGTH = ",i4,1x,a)') j,cnt,temp
    ELSE 
      WRITE(UNIT=num,FMT='(6x,"SEQUENCE ",i2,": DNA     LENGTH = ",i4,1x,a)') j,cnt,temp
    END IF

    WRITE(UNIT=num,FMT="(1x,a)") bar64

    nlines=cnt/60                    ! 60 aa/line

! Write the sequence with 60 nt per line

    DO i=1,nlines
      ne=i*60
      ns=ne-59
      WRITE(UNIT=num,FMT="(i4,1x,a60)") ns,INITseq(ns+start-1:ne+start-1)
    END DO

! Now write the last line of the sequence (less than 60 nt)

    ns=nlines*60+1                  ! find the current nt number
    ne=cnt                      ! find the total length of sequence
    WRITE(UNIT=num,FMT="(i4,1x,a)") ns,INITseq(ns+start-1:ne+start-1)
    WRITE(UNIT=num,FMT="(1x,a)") bar64
  END DO

END SUBROUTINE Print_Seq_Log
SUBROUTINE Print_SimpleDNA(num)
!
! Prints a the current DNA sequence in simple format to the console
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,ns,ne,nlines,num

  IF (TEST0) PRINT *,"Print_SimpleDNA" !TEST0

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT='(16x,"INITIAL DNA SEQUENCE")')
  WRITE(UNIT=num,FMT="(1x,a64)") bar64

  nlines=DNAlen/60
  DO i=1,nlines
    ne=i*60
    ns=ne-59
    WRITE(UNIT=num,FMT="(i4,1x,a60)") ns,CurrDNA%DNAseq(ns:ne)
  END DO
  ns=nlines*60+1                  ! find the current nt number
  ne=DNAlen                       ! find the total length of sequence
  WRITE(UNIT=num,FMT="(i4,1x,a/)") ns,CurrDNA%DNAseq(ns:ne)

END SUBROUTINE Print_SimpleDNA
SUBROUTINE Print_SimplePROT(num)
!
! Prints a the current protein sequence in simple format to the console
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,ns,ne,nlines,num

  IF (TEST0) PRINT *,'Print_SimplePROT'

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT='(16x,"INITIAL PROTEIN SEQUENCE")')
  WRITE(UNIT=num,FMT="(1x,a64)") bar64

  nlines=PROTlen/60
  DO i=1,nlines
    ne=i*60
    ns=ne-59
    WRITE(UNIT=num,FMT="(i4,1x,a60)") ns,PROTseq(ns:ne)
  END DO
  ns=nlines*60+1                  ! find the current nt number
  ne=PROTlen                       ! find the total length of sequence
  WRITE(UNIT=num,FMT="(i4,1x,a/)") ns,PROTseq(ns:ne)

END SUBROUTINE Print_SimplePROT
SUBROUTINE Print_TranslatedDNA(num)
!
! Prints DNA with protein sequence
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,ns,ne,nlines,num
  CHARACTER(LEN=80) :: line

  IF (TEST0) PRINT *,'Print_TranslatedDNA'

  WRITE(UNIT=num,FMT="(' ')")
  WRITE(UNIT=num,FMT='(19x,"INITIAL TRANSLATED SEQUENCE")')
  WRITE(UNIT=num,FMT="(1x,a64)") bar64

  nlines=DNAlen/60
  DO i=1,nlines

! dna sequence

    line=''
    ne=i*60
    ns=ne-59
    WRITE(line,FMT="(i4)") ns     ! first nt number of line
    DO j=0,59
      line(j+6:j+6)=CurrDNA%DNAseq(j+ns:j+ns)
    END DO
    WRITE(UNIT=num,FMT="(a80)") line

! complement

    CALL ComplStr(line(6:65))
    WRITE(UNIT=num,FMT="(a65)") line(6:65)

! protein sequence

    line=''
    DO j=0,59
      IF (nt2aa(j+ns).ne.0) line(j+6:j+6)=AAT(nt2aa(j+ns))%AA1
    END DO
    WRITE(UNIT=num,FMT="(a80)") line

! blank line

    WRITE(UNIT=num,FMT="('')") 
  END DO

! last line

  ns=nlines*60+1                  ! find the current nt number
  ne=DNAlen                       ! find the total length of sequence
  WRITE(line,FMT="(i4)") ns     ! first nt number of line
  DO j=0,59
    line(j+6:j+6)=CurrDNA%DNAseq(j+ns:j+ns)
  END DO
  WRITE(UNIT=num,FMT="(a80)") line              ! dna sequence
  CALL ComplStr(line(6:65))
  WRITE(UNIT=num,FMT="(a65)") line(6:65)        ! complement

! protein sequence

  line=''
  DO j=0,59
    IF (nt2aa(j+ns).ne.0) line(j+6:j+6)=AAT(nt2aa(j+ns))%AA1
  END DO
  WRITE(UNIT=num,FMT="(a80)") line

! blank line

  WRITE(UNIT=num,FMT="('')") 

END SUBROUTINE Print_TranslatedDNA
