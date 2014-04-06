SUBROUTINE Add_DNA(start,finish,reverse,gf)
!
! Adds a block of nucleotides to the growing sequence
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,start,finish,first
  CHARACTER(LEN=1) :: nt
  INTEGER,EXTERNAL :: NT2Int
  LOGICAL :: reverse,done,gf
  CHARACTER(LEN=9999) :: tempDNAseq

  done=.FALSE.
  IF (TEST0) PRINT *,"Add_DNA" !TEST0

  NumberOfSeq=NumberOfSeq+1

! degenerate sequences are allowed for nucleotides... to a point

  DO i=start,finish
    DO j=1,LEN_TRIM(InputArrayUC(i))
      nt=InputArrayUC(i)(j:j)
      IF (nt.eq.'A'.or.nt.eq.'C'.or.nt.eq.'G'.or.nt.eq.'T'.or.nt.eq.'N'.or. &
          nt.eq.'M'.or.nt.eq.'R'.or.nt.eq.'W'.or.nt.eq.'S'.or.nt.eq.'Y'.or. &
          nt.eq.'K'.or.nt.eq.'B'.or.nt.eq.'D'.or.nt.eq.'H'.or.nt.eq.'V') THEN
        DNAlen=DNAlen+1
        INITlen=INITlen+1
        INIT2Seq(INITlen)=NumberOfSeq
        IF (DNAlen.gt.MaxDNAlen) CALL Stop_Program("Too many nucleotides.")
        CurrDNA%DNAseq(DNAlen:DNAlen)=nt
        INITseq(INITlen:INITlen)=nt
        ORIGDNAseq(DNAlen:DNAlen)=nt
        CurrDNA%NUMseq(DNAlen)=NT2Int(nt) ! if the nt is degenerate, this is undefined 
        IF (.not.done) THEN
          first = DNAlen
          done = .TRUE.
        END IF
      END IF
      IF (nt.eq.'M'.or.nt.eq.'R'.or.nt.eq.'W'.or.nt.eq.'S'.or.nt.eq.'Y'.or. & 
          nt.eq.'K'.or.nt.eq.'B'.or.nt.eq.'D'.or.nt.eq.'H'.or.nt.eq.'V'.or. & 
          nt.eq.'N') THEN
        IF (.not.gf) CALL Stop_Program("Degenerate sequences must be gapfixed.")
        NumDegPos=NumDegPos+1
        DegPos(NumDegPos)=DNAlen
        CurrDNA%Degen(DNAlen) = .TRUE.

! determine degenerate index and assign it
        DO k=1,11
          IF (DegenSeq(k)%DegNT.eq.nt) THEN
            CurrDNA%DegenNum(DNAlen) = k
          END IF
        END DO

      ELSE
        CurrDNA%Degen(DNAlen) = .FALSE.
      END IF
    END DO
  END DO

! reverse the nucleotide sequence if desired
  IF (reverse) THEN
    SeqReverse(NumberOfSeq)=.TRUE.
    CALL RevComplStr(CurrDNA%DNAseq(first:DNAlen))
    DO i=first,DNAlen
      CurrDNA%NUMseq(i)=NT2Int(CurrDNA%DNAseq(i:i))
    END DO
  END IF

! assign the GapFix array if the positions need to be in a gap
  IF (gf) THEN
    SeqGapFix(NumberOfSeq)=.TRUE.
    DO i=first,DNAlen
      CurrDNA%GapFixPos(i)=.TRUE.
    END DO
  END IF

END SUBROUTINE Add_DNA
SUBROUTINE Add_Previous

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=256) :: a1,a2,a3,a4,a5,a6 

  INTEGER,EXTERNAL :: StrToInt
  REAL,EXTERNAL :: StrToReal
  LOGICAL,EXTERNAL :: Read_Old_Logfile

  INTEGER :: ierr,i,j,k,i1,i2 

  IF (TEST0) PRINT *,"Add_Previous" !TEST0

! PREVIOUS
! MutantRun
! PrevTrial

  prev: DO i=1,InputArrayNum
    a1=""
    a2=""
    a3=""
    IF (INDEX(InputArrayUC(i),'PREVIOUS').eq.1) THEN
      READ(InputArray(i),*,IOSTAT=ierr) a1,a2,a3  ! are there three lines?
      IF (ierr.eq.0) THEN
        IF ((LEN_TRIM(a3)).ne.0) oldlogfile=a3
      ELSE
        READ(InputArrayUC(i),*,IOSTAT=ierr) a1,a2
        IF (ierr.ne.0) EXIT prev
      END IF
      PrevTrial=StrToInt(a2)
      IF (PrevTrial.gt.0) THEN
        MutantRun=.TRUE.

! try to read the oldlogfile -- if successful, replace all the relevant info

        IF (Read_Old_Logfile(PrevTrial)) THEN
          FinalScore(1)%Oligo=OligoLen
          FinalScore(1)%MeltT=MeltTemp
          TotalNumberOfSolutions=1
          WRITE(a4,FMT='(i3)') PrevTrial
          a4=TRIM(a4)
          a4=ADJUSTL(a4)
          i2=LEN_TRIM(jobname)
          j=LEN_TRIM(OLDjobname)
          k=LEN_TRIM(a4)
          jobname=jobname(1:i2)//' (using trial '//a4(1:k)//' parameters from previous job '//OLDjobname(1:j)//')'
        END IF
      ELSE
        PrevTrial=0
      END IF
    EXIT prev
    END IF
  END DO prev

END SUBROUTINE Add_Previous
SUBROUTINE Add_Protein(start,finish,reverse,gf)
!
! Adds residue to PROTseq, translates protein, fills in growing arrays
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,n,x,start,finish,z,first,last,tempDNAlen,tempPROTlen
  INTEGER :: pos1,pos2,pos3
  REAL :: rand
  LOGICAL :: no_codons
  LOGICAL :: done
  LOGICAL :: reverse,gf
  INTEGER,EXTERNAL :: NT2Int
  CHARACTER(LEN=3) :: tempCodonSeq

  done=.FALSE.
  tempDNAlen=0
  tempPROTlen=0
  
  IF (TEST0) PRINT *,"Add_Protein" !TEST0

! first get the basic sequence and amino acid type

  DO i=start,finish
    DO z=1,LEN_TRIM(InputArrayUC(i))
      chk1: DO j=1,21
        IF (InputArrayUC(i)(z:z).eq.AAT(j)%AA1) THEN

! Add the residue to the protein sequence

          PROTlen=PROTlen+1
          tempPROTlen=tempPROTlen+1
          tempDNAlen=tempDNAlen+3
          IF (PROTlen.gt.MaxPROTlen) CALL Stop_Program("Too many protein residues.")
          PROTseq(PROTlen:PROTlen)=AAT(j)%AA1
          prot2aa(PROTlen)=j

! Build up the initial sequence for logfile output

          INITlen=INITlen+1
          INITseq(INITlen:INITlen)=AAT(j)%AA1

! Increment NumberOfChains and prot2chain array

          IF (.not.done) THEN
            NumberOfChains=NumberOfChains+1
            NumberOfSeq=NumberOfSeq+1
            first = PROTlen
            IF (reverse) ChainReverse(NumberOfChains)=.TRUE.
            IF (reverse) SeqReverse(NumberOfSeq)=.TRUE.
            IF (gf) ChainGapFix(NumberOfChains)=.TRUE.
            IF (gf) SeqGapFix(NumberOfSeq)=.TRUE.
            SeqIsProt(NumberOfSeq)=.TRUE.
          END IF
          prot2chain(PROTlen)=NumberOfChains
          INIT2Seq(INITlen)=NumberOfSeq
  
! Next residue, please

          done=.TRUE. 

          EXIT chk1
        END IF
      END DO chk1
    END DO
  END DO

  last=PROTlen
  j=0

! now create the links between DNA sequence and PROT sequence

  DO i=first,last

    j=j+3

! Either choose the codon with the highest frequency

    k = 1 

! or choose the codon randomly unless the codon is not allowed

    IF (CodonRandom) THEN
      CALL RANDOM_NUMBER(rand)
      k=(INT(rand*(AAT(prot2aa(i))%NumOfActiveCodons)))+1
    END IF

    tempCodonSeq=CFT(AAT(prot2aa(i))%Codon(k))%Seq

! Increment mutPROT2prot array
      
    IF (AAT(prot2aa(i))%NumOfActiveCodons.gt.1) THEN
      mutPROTnum = mutPROTnum+1
      mutPROT2prot(mutPROTnum) = i
    END IF

! Fill prot2cod arrays

    CurrDNA%prot2cod(i) = AAT(prot2aa(i))%Codon(k)

! If reverse is true, put protein in backward

    IF (.not.reverse) THEN
      pos1=DNAlen+j-2
      pos2=DNAlen+j-1
      pos3=DNAlen+j      
    ELSE 
      pos1=DNAlen+tempDNAlen-j+1
      pos2=DNAlen+tempDNAlen-j+2
      pos3=DNAlen+tempDNAlen-j+3
      CALL RevComplStr(tempCodonSeq)
    END IF

! Fill prot2nt, nt2aa, nt2prot, and prot2nt arrays

    nt2aa(pos2) = prot2aa(i)
    prot2nt(i) = pos2
    nt2prot(pos1) = i
    nt2prot(pos2) = i
    nt2prot(pos3) = i

! Insert codon into the DNA sequence.

    CurrDNA%DNAseq(pos1:pos3)=tempCodonSeq

! Fill nt2codon arrays

    CurrDNA%nt2cod(pos2) = AAT(prot2aa(i))%Codon(k)

! Fill the numerical sequence array

    CurrDNA%NUMseq(pos1)=NT2Int(CurrDNA%DNAseq(pos1:pos1))
    CurrDNA%NUMseq(pos2)=NT2Int(CurrDNA%DNAseq(pos2:pos2))
    CurrDNA%NUMseq(pos3)=NT2Int(CurrDNA%DNAseq(pos3:pos3))

! Fill the GapFixPos array

    IF (gf) THEN
      CurrDNA%GapFixPos(pos1) = .TRUE.
      CurrDNA%GapFixPos(pos2) = .TRUE.
      CurrDNA%GapFixPos(pos3) = .TRUE.
    END IF

  END DO

  DNAlen=DNAlen+tempDNAlen

END SUBROUTINE Add_Protein
SUBROUTINE Default_Param
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k
  CHARACTER(LEN=63) :: a,b
  CHARACTER(LEN=192) :: c
  INTEGER,EXTERNAL :: CurrentTimeSeconds
  INTEGER :: values(8)
  CHARACTER(LEN=8) :: date
  CHARACTER(LEN=10) :: time
  CHARACTER(LEN=5) :: zone
  INTEGER,EXTERNAL :: NT2Int

  IF (TEST0) PRINT *,"Default_Param" !TEST0

  a='ACDEFGHIKLMNPQRSTVWXY                                          '
  b='AlaCysAspGluPheGlyHisIleLysLeuMetAsnProGlnArgSerThrValTrpEndTyr'
  c='TTTTCTTATTGTTTCTCCTACTGCTTATCATAATGATTGTCGTAGTGGCTTCCTCATCGTCTCC&
    &CCCACCGCCTACCACAACGACTGCCGCAGCGGATTACTAATAGTATCACCAACAGCATAACAAA&
    &AAGAATGACGAAGAGGGTTGCTGATGGTGTCGCCGACGGCGTAGCAGAAGGAGTGGCGGAGGGG'

! Initialize the arrays

  DO i=1,9999
    nt2aa(i) = 0
    CurrDNA%nt2cod(i) = 0
    CurrDNA%GapFixPos(i) = .FALSE.
  END DO

  DO i=1,21
    j=((i-1)*3)+1
    AAT(i)%AA1=a(i:i)                  ! Assign One Letter Aa values
    AAT(i)%AA3=b(j:j+2)              ! Assign Three Letter Aa values
    AAT(i)%NumOfCodons=0
    AAT(i)%NumOfActiveCodons=0
    AAT(i)%NumberSum=0
    DO k=1,10
      AAT(i)%Codon(k)=0
      AAT(i)%Freq(k)=0
    END DO
  END DO

  DO i=1,64
    j=((i-1)*3)+1
    CFT(i)%Seq=c(j:j+2)
    CFT(i)%SeqRC(1:3)=c(j:j+2)
    CALL RevComplStr(CFT(i)%SeqRC(1:3))
    CFT(i)%Check=.TRUE.                    ! all are true by default
    CFT(i)%AA1=''
    CFT(i)%AA3=''
    CFT(i)%Freq=0
    CFT(i)%Number=0
    DO k=1,3
      CFT(i)%num(k)=NT2Int(CFT(i)%Seq(k:k))
      CFT(i)%numRC(k)=NT2Int(CFT(i)%SeqRC(k:k))
    END DO
  END DO

  DO i=1,99
    ChainReverse(i)=.FALSE.
    ChainGapFix(i)=.FALSE.
    SeqReverse(i)=.FALSE.
    SeqGapFix(i)=.FALSE.
    SeqIsProt(i)=.FALSE.
  END DO

! Degenerate sequences

  DegenSeq(1)%DegNT = 'K'
  DegenSeq(1)%NumOfNT = 2
  DegenSeq(1)%Seq(1) = 'G'
  DegenSeq(1)%NumSeq(1) = NT2Int(DegenSeq(1)%Seq(1))
  DegenSeq(1)%Seq(2) = 'T'
  DegenSeq(1)%NumSeq(2) = NT2Int(DegenSeq(1)%Seq(2))

  DegenSeq(2)%DegNT = 'M'
  DegenSeq(2)%NumOfNT = 2
  DegenSeq(2)%Seq(1) = 'A'
  DegenSeq(2)%NumSeq(1) = NT2Int(DegenSeq(2)%Seq(1))
  DegenSeq(2)%Seq(2) = 'C'
  DegenSeq(2)%NumSeq(2) = NT2Int(DegenSeq(2)%Seq(2))

  DegenSeq(3)%DegNT = 'R'
  DegenSeq(3)%NumOfNT = 2
  DegenSeq(3)%Seq(1) = 'A'
  DegenSeq(3)%NumSeq(1) = NT2Int(DegenSeq(3)%Seq(1))
  DegenSeq(3)%Seq(2) = 'G'
  DegenSeq(3)%NumSeq(2) = NT2Int(DegenSeq(3)%Seq(2))

  DegenSeq(4)%DegNT = 'S'
  DegenSeq(4)%NumOfNT = 2
  DegenSeq(4)%Seq(1) = 'C'
  DegenSeq(4)%NumSeq(1) = NT2Int(DegenSeq(4)%Seq(1))
  DegenSeq(4)%Seq(2) = 'G'
  DegenSeq(4)%NumSeq(2) = NT2Int(DegenSeq(4)%Seq(2))

  DegenSeq(5)%DegNT = 'W'
  DegenSeq(5)%NumOfNT = 2
  DegenSeq(5)%Seq(1) = 'A'
  DegenSeq(5)%NumSeq(1) = NT2Int(DegenSeq(5)%Seq(1))
  DegenSeq(5)%Seq(2) = 'T'
  DegenSeq(5)%NumSeq(2) = NT2Int(DegenSeq(5)%Seq(2))

  DegenSeq(6)%DegNT = 'Y'
  DegenSeq(6)%NumOfNT = 2
  DegenSeq(6)%Seq(1) = 'C'
  DegenSeq(6)%NumSeq(1) = NT2Int(DegenSeq(6)%Seq(1))
  DegenSeq(6)%Seq(2) = 'T'
  DegenSeq(6)%NumSeq(2) = NT2Int(DegenSeq(6)%Seq(2))

  DegenSeq(7)%DegNT = 'B'
  DegenSeq(7)%NumOfNT = 3
  DegenSeq(7)%Seq(1) = 'C'
  DegenSeq(7)%NumSeq(1) = NT2Int(DegenSeq(7)%Seq(1))
  DegenSeq(7)%Seq(2) = 'G'
  DegenSeq(7)%NumSeq(2) = NT2Int(DegenSeq(7)%Seq(2))
  DegenSeq(7)%Seq(3) = 'T'
  DegenSeq(7)%NumSeq(3) = NT2Int(DegenSeq(7)%Seq(3))

  DegenSeq(8)%DegNT = 'D'
  DegenSeq(8)%NumOfNT = 3
  DegenSeq(8)%Seq(1) = 'A'
  DegenSeq(8)%NumSeq(1) = NT2Int(DegenSeq(8)%Seq(1))
  DegenSeq(8)%Seq(2) = 'C'
  DegenSeq(8)%NumSeq(2) = NT2Int(DegenSeq(8)%Seq(2))
  DegenSeq(8)%Seq(3) = 'T'
  DegenSeq(8)%NumSeq(3) = NT2Int(DegenSeq(8)%Seq(3))

  DegenSeq(9)%DegNT = 'H'
  DegenSeq(9)%NumOfNT = 3
  DegenSeq(9)%Seq(1) = 'A'
  DegenSeq(9)%NumSeq(1) = NT2Int(DegenSeq(9)%Seq(1))
  DegenSeq(9)%Seq(2) = 'C'
  DegenSeq(9)%NumSeq(2) = NT2Int(DegenSeq(9)%Seq(2))
  DegenSeq(9)%Seq(3) = 'T'
  DegenSeq(9)%NumSeq(3) = NT2Int(DegenSeq(9)%Seq(3))

  DegenSeq(10)%DegNT = 'V'
  DegenSeq(10)%NumOfNT = 3
  DegenSeq(10)%Seq(1) = 'A'
  DegenSeq(10)%NumSeq(1) = NT2Int(DegenSeq(10)%Seq(1))
  DegenSeq(10)%Seq(2) = 'C'
  DegenSeq(10)%NumSeq(2) = NT2Int(DegenSeq(10)%Seq(2))
  DegenSeq(10)%Seq(3) = 'G'
  DegenSeq(10)%NumSeq(3) = NT2Int(DegenSeq(10)%Seq(3))

  DegenSeq(11)%DegNT = 'N'
  DegenSeq(11)%NumOfNT = 4
  DegenSeq(11)%Seq(1) = 'A'
  DegenSeq(11)%NumSeq(1) = NT2Int(DegenSeq(11)%Seq(1))
  DegenSeq(11)%Seq(2) = 'C'
  DegenSeq(11)%NumSeq(2) = NT2Int(DegenSeq(11)%Seq(2))
  DegenSeq(11)%Seq(3) = 'G'
  DegenSeq(11)%NumSeq(3) = NT2Int(DegenSeq(11)%Seq(3))
  DegenSeq(11)%Seq(4) = 'T'
  DegenSeq(11)%NumSeq(4) = NT2Int(DegenSeq(11)%Seq(4))

END SUBROUTINE Default_Param
SUBROUTINE Fix_Codons
!
! Determines if 64 codons are present, sorts frequencies, and assigns check
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,n
  REAL :: x,y

  IF (TEST0) PRINT *,"Fix_Codons" !TEST0

! Stop if not enough codons are read in

  j=0
  DO i=1,64
    IF (CFT(i)%Freq.ne.0) j=j+1
  END DO
  IF (j.ne.64) CALL Stop_Program("There are too few codons.")

! Sort the frequencies from highest to lowest for each AA
! (i.e., AAT(i)%Freq(1) = highest, AAT(i)%Freq(..NumOfCodons) = lowest

! This is a change --
! To sort from lowest to highest, switch lt to gt in comparison  7/27/05

  DO i=1,21
    IF (AAT(i)%NumOfCodons.gt.1) THEN
      j=AAT(i)%NumOfCodons-1
      DO n=j,1,-1
        DO k=1,n
          IF (AAT(i)%Freq(k).lt.AAT(i)%Freq(k+1)) THEN
            CALL RealSwap(AAT(i)%Freq(k),AAT(i)%Freq(k+1))
            CALL IntSwap(AAT(i)%Codon(k),AAT(i)%Codon(k+1))
          END IF
        END DO
      END DO
    END IF
  END DO
   
! Assign CFT%Check values according to SeqOptimTol

  DO i=1,21
    DO j=1,AAT(i)%NumOfCodons

! If the SeqOptimToler is higher than the particular codon frequency for that
! amino acid, deactivate that codon

! Unless SeqOptimToler equals 100, which allows scoring and all codons to
! be used.

      x=AAT(i)%Freq(j)*100
      y=REAL(SeqOptimToler)
      IF (SeqOptimToler.lt.100) THEN
        IF (AAT(i)%Freq(j)*100.lt.SeqOptimToler) THEN
          CFT(AAT(i)%Codon(j))%Check=.FALSE.
        END IF
      END IF
    END DO

! However, keep the two highest frequency codons available for generating DNA

    CFT(AAT(i)%Codon(1))%Check=.TRUE.
    IF (AAT(i)%NumofCodons.gt.1.and.(.not.CodonStrict)) THEN
      CFT(AAT(i)%Codon(2))%Check=.TRUE.
    END IF

  END DO

! Assign NumOfActiveCodons

  DO i=1,21
    DO j=1,AAT(i)%NumOfCodons
      IF (CFT(AAT(i)%Codon(j))%Check) THEN
        AAT(i)%NumOfActiveCodons=AAT(i)%NumOfActiveCodons+1
      END IF
    END DO
  END DO

END SUBROUTINE Fix_Codons
SUBROUTINE Read_Codon_Line(AA3,codon,y,x)
!
! This subroutine reads in a single line of codon frequency file (GCG format).

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k
  REAL :: freq,number
  CHARACTER(LEN=5) :: AA3
  CHARACTER(LEN=5) :: codon
  CHARACTER(LEN=5) :: x 
  CHARACTER(LEN=12) :: y
  REAL,EXTERNAL :: StrToReal

  IF (TEST1) PRINT *,"Read_Codon_Line"

  IF (AA3.ne.''.and.codon.ne.''.and.x.ne.'') THEN

! fix the input

    AA3=ADJUSTL(AA3)
    codon=ADJUSTL(codon)
    freq=StrToReal(x)
    number=StrToReal(y)

    IF (AA3(1:3).eq."Stp") AA3(1:3) = "End"
    IF (AA3(1:3).eq."Sto") AA3(1:3) = "End"
    IF (AA3(1:3).eq."STO") AA3(1:3) = "End"

! Load the temporary variables into program input variables

    loop2: DO j=1,64
      IF (codon(1:3).eq.CFT(j)%Seq) THEN
        loop3: DO k=1,21
          IF (AA3(1:3).eq.AAT(k)%AA3) THEN
            AAT(k)%NumOfCodons=AAT(k)%NumOfCodons+1
            AAT(k)%Codon(AAT(k)%NumOfCodons)=j
            AAT(k)%Freq(AAT(k)%NumOfCodons)=freq
            AAT(k)%NumberSum=AAT(k)%NumberSum+number
            CFT(j)%AA3 = AA3(1:3)
            CFT(j)%AA1 = AAT(k)%AA1
            CFT(j)%Freq = freq
            CFT(j)%Number = number
            EXIT loop3
          END IF
        END DO loop3
        EXIT loop2
      END IF
    END DO loop2

  END IF


END SUBROUTINE Read_Codon_Line
SUBROUTINE Read_Codons

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=256) :: a1,a2,a3,a4,a5,a6 ! throw away strings

  INTEGER,EXTERNAL :: StrToInt
  REAL,EXTERNAL :: StrToReal

  INTEGER :: ierr,length,i,j,k,n,start,finish
  INTEGER :: i1,i2,i3,i4,i5,i6  ! throw away integers

  IF (TEST0) PRINT *,"Read_Codons" !TEST0

  main: DO k=1,InputArrayNum

! read an individual line, trim it, and change to upper case

    a1=""
    a2=""
    a3=""
    a4=""
    a5=""

! CODON frequency table

    a1=InputArrayUC(k)
    IF (INDEX(a1,'CODON').eq.1) THEN
      IF ((INDEX(a1,'ECOLI2')).gt.0) THEN
        CALL Read_Organism_CFT('ecoli2')
      ELSE IF ((INDEX(a1,'COLI')).gt.0.or.(INDEX(a1,'ESCHER')).gt.0) THEN
        CALL Read_Organism_CFT('ecoli ')
      ELSE IF ((INDEX(a1,'ELEGAN')).gt.0.or.(INDEX(a1,'CAENOR')).gt.0) THEN
        CALL Read_Organism_CFT('celega')
      ELSE IF ((INDEX(a1,'MELANO')).gt.0.or.(INDEX(a1,'DROSOP')).gt.0) THEN
        CALL Read_Organism_CFT('drosop')
      ELSE IF ((INDEX(a1,'SAPIEN')).gt.0.or.(INDEX(a1,'HOMO S')).gt.0) THEN
        CALL Read_Organism_CFT('homosa')
      ELSE IF ((INDEX(a1,'MUSCUL')).gt.0.or.(INDEX(a1,'MUS MU')).gt.0) THEN
        CALL Read_Organism_CFT('musmus')
      ELSE IF ((INDEX(a1,'PASTOR')).gt.0.or.(INDEX(a1,'PICHIA')).gt.0) THEN
        CALL Read_Organism_CFT('pastor')
      ELSE IF ((INDEX(a1,'NORVEG')).gt.0.or.(INDEX(a1,'RATTUS')).gt.0) THEN
        CALL Read_Organism_CFT('rattus')
      ELSE IF ((INDEX(a1,'CEREVE')).gt.0.or.(INDEX(a1,'SACCHA')).gt.0) THEN
        CALL Read_Organism_CFT('saccho')
      ELSE IF ((INDEX(a1,'LAEVIS')).gt.0.or.(INDEX(a1,'XENOPU')).gt.0) THEN
        CALL Read_Organism_CFT('xenopu')
      ELSE
        Organism="USER INPUT"

! Find end of codon block

        start=k+1
        seek: DO j=k,InputArrayNum
          IF (INDEX(InputArrayUC(j),'//').gt.0) THEN
            finish=j-1
            EXIT seek ! end of codon block
          END IF
        END DO seek
     
        codo: DO j=start,finish
          READ(InputArray(j),*,IOSTAT=ierr) a1,a2,a3,a4,a5
          IF (ierr.ne.0) CYCLE codo
          CALL Read_Codon_Line(a1,a2,a3,a5)
        END DO codo

      END IF
    EXIT main
    END IF

  END DO main

! Make sure frequencies are set

  k=0
  DO i=1,64
    IF (CFT(i)%Freq.eq.0) THEN
      k=k+1
    END IF
  END DO

  IF (k.eq.64) THEN
    DO i=1,21
      DO j=1,AAT(i)%NumOfCodons
        AAT(i)%Freq(j)=CFT(AAT(i)%Codon(j))%Number/AAT(i)%NumberSum
        CFT(AAT(i)%Codon(j))%Freq=CFT(AAT(i)%Codon(j))%Number/AAT(i)%NumberSum
      END DO
    END DO
  END IF
      

! Make sure the input makes sense

  k=0
  DO i=1,21
    k=k+AAT(i)%NumOfCodons
  END DO

! Default CFT is E. coli Class II

  IF (k.ne.64) THEN
    DO i=1,21
      AAT(i)%NumOfCodons=0
      AAT(i)%NumOfActiveCodons=0
    END DO
    CALL Read_Organism_CFT('ecoli2')
  END IF

  CALL Fix_Codons

END SUBROUTINE Read_Codons
SUBROUTINE Read_Input

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

! Here are the things needed for input:
!
! OligoLenHi      user input oligo size (upper limit)
! OligoLenLo      user input oligo size (lower limit)
! MeltTempHi      Ideal melting temperature (upper limit)
! MeltTempLo      Ideal melting temperature (lower limit)
! MeltTol         Tolerance for melting temperature deviation
! SeqOptimToler   Lowest allowed codon frequency
! NumberOfSolutions  
! OligoConcTemp   in nM
! SodiumConc      in mM
! MgConc          in mM
! TBIOSwitch      off or on
! RepLen          determines the size of repeats to minimize
! MPLn            length of  misprimes
! MaxNonId        maximum number of non-identical nts in misprime
! MPTip           identical tip of the misprime, in nts
! PrevTrial       number of previous run

  CHARACTER(LEN=256) :: line
  CHARACTER(LEN=256) :: a1,a2,a3,a4,a5,a6 ! throw away strings
  CHARACTER(LEN=3) :: temp

  INTEGER,EXTERNAL :: StrToInt
  REAL,EXTERNAL :: StrToReal
  LOGICAL,EXTERNAL :: Read_Old_Logfile
  LOGICAL :: reverse
  LOGICAL :: gf
  LOGICAL :: skipdquote,skipsquote

  INTEGER :: ierr,length,i,j,k,n,start,finish
  INTEGER :: i1,i2,i3,i4,i5,i6  ! throw away integers

  IF (TEST0) PRINT *,"Read_Input" !TEST0

  OPEN (UNIT=inputnum,FILE=inputfile,STATUS='old',IOSTAT=ierr)
  IF (ierr.ne.0) THEN
    WRITE(UNIT=console,FMT="('')")
    WRITE(UNIT=console,FMT="('Where is the inputfile ""',a,'""?')") inputfile(1:LEN_TRIM(inputfile))
    PRINT *,""
    PRINT *,"Usage: dnaworks [ inputfile ] [ -help ]"
    PRINT *,""
    STOP
  END IF

! Read the input file into the InputArray array

  i1=0
  reading: DO 
    READ(UNIT=inputnum,FMT='(a)',IOSTAT=ierr) line
    IF (ierr.ne.0) EXIT reading     ! dump out if end of file reached
    i1=i1+1

! Take out comments

    IF ((INDEX(line,'!')).gt.0) line(INDEX(line,'!'):256)=''
    IF ((INDEX(line,'#')).gt.0) line(INDEX(line,'#'):256)=''

    InputArray(i1) = line

! Take out quoted strings and create uppercase array as well

    skipsquote=.FALSE.
    skipdquote=.FALSE.
    DO i=1,LEN_TRIM(line)
      IF (skipsquote) THEN
        IF (line(i:i).eq."'") THEN
          line(i:i)=' '
          skipsquote=.FALSE.
        ELSE
          line(i:i)=' '
        END IF
      ELSE
        IF (line(i:i).eq."'") THEN
          line(i:i)=' '
          skipsquote=.TRUE. 
        END IF
      END IF
      IF (skipdquote) THEN
        IF (line(i:i).eq.'"') THEN
          line(i:i)=' '
          skipdquote=.FALSE.
        ELSE
          line(i:i)=' '
        END IF
      ELSE
        IF (line(i:i).eq.'"') THEN
          line(i:i)=' '
          skipdquote=.TRUE. 
        END IF
      END IF
    END DO

    InputArrayUC(i1) = line

    CALL ToUpperCase(InputArrayUC(i1))

  END DO reading
  CLOSE (UNIT=inputnum)
  InputArrayNum=i1

! It is very important to read these in order... 

  CALL Read_Parameters
  CALL Read_Codons
  CALL Read_Patterns
  
! step through the input

  main: DO i=1,InputArrayNum

! NUCLEOTIDE

    IF (INDEX(InputArrayUC(i),'NUCLEOTIDE').eq.1) THEN
      reverse=.FALSE.
      gf=.FALSE.
      IF (INDEX(InputArrayUC(i),'REVERSE').gt.0) reverse=.TRUE.
      IF (INDEX(InputArrayUC(i),'GAPFIX').gt.0) gf=.TRUE.
      IF (INDEX(InputArrayUC(i),'GAPFIX').gt.0) GapFix=.TRUE.
      start=i+1
      seek1: DO j=i,InputArrayNum
        IF (INDEX(InputArray(j),'//').gt.0) THEN
          finish=j
          EXIT seek1
        END IF
      END DO seek1 
      CALL Add_DNA(start,finish,reverse,gf)
    END IF

! PROTEIN

    IF (INDEX(InputArrayUC(i),'PROTEIN').eq.1) THEN
      reverse=.FALSE.
      gf=.FALSE.
      IF (INDEX(InputArrayUC(i),'REVERSE').gt.0) reverse=.TRUE.
      IF (INDEX(InputArrayUC(i),'GAPFIX').gt.0) gf=.TRUE.
      IF (INDEX(InputArrayUC(i),'GAPFIX').gt.0) GapFix=.TRUE.
      start=i+1
      seek2: DO j=i,InputArrayNum
        IF (INDEX(InputArrayUC(j),'//').gt.0) THEN
          finish=j
          EXIT seek2
        END IF
      END DO seek2
      CALL Add_Protein(start,finish,reverse,gf)
    END IF

  END DO main

! After reading in everything else, find if prev trial is set

  CALL Add_Previous

! Don't forget to create OligoCorr and SaltCorr!

  CALL TmCorrect

END SUBROUTINE Read_Input
SUBROUTINE Read_Old_Codons

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: ierr,i,j,k,m,z
  CHARACTER(LEN=100) :: line
  INTEGER,EXTERNAL :: StrToInt
  REAL,EXTERNAL :: StrToReal
  TYPE(KnownCodon) :: tempCFT(64)  ! Codon Frequency Table
  TYPE(KnownAA) :: tempAAT(21)     ! Amino Acid Table
  CHARACTER(LEN=10) :: a1(4),a2(4)
  REAL :: r1(4),r2(4)
  INTEGER :: i1(4),i2(4)
  LOGICAL :: oldCodons=.FALSE.     ! are old codons present?

  IF (TEST0) PRINT *,'Read_Old_Codons'  

  OPEN (UNIT=oldlognum,FILE=oldlogfile,STATUS='old',IOSTAT=ierr)

! -- CODON FREQUENCY TABLE --

  main: DO
    READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
    IF (ierr.ne.0) EXIT main
    IF (INDEX(line,'CODON FREQUENCY TABLE').gt.0) THEN
!      PRINT *,'codon input'
      oldCodons=.TRUE.
      Organism=line((LogfileOffset+36):(LogfileOffset+65))
      m=0
      codon: DO z=1,22
        READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
        IF (ierr.ne.0) EXIT main
        READ(line,*,IOSTAT=ierr) a1(1),a2(1),r1(1),a1(2),a2(2),r1(2),a1(3),a2(3),r1(3),a1(4),a2(4),r1(4)
        IF (ierr.ne.0) CYCLE codon
        IF (m.eq.64) EXIT codon 
        DO i=1,4
          m=m+1
          tempCFT(m)%Seq=a1(i)(1:3)
          tempCFT(m)%AA1=a2(i)(1:1)
          IF (r1(i).le.0) r1(i)=0.001
          tempCFT(m)%Freq=r1(i)
!          PRINT *,m,tempCFT(m)%Seq,tempCFT(m)%Freq
        END DO
      END DO codon
    END IF

! pull out other AAT and CFT information

! -- ACTIVE CODONS --

    IF (INDEX(line,'ACTIVE CODONS').gt.0) THEN
!      PRINT *,'active codons'
      m=0
      amino: DO i=1,31
        READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
        IF (ierr.ne.0) EXIT main
        READ(line,*,IOSTAT=ierr) a1(1),a2(1),i1(1),i2(1),r1(1),r2(1)
        IF (ierr.ne.0) CYCLE amino
        m=m+1
        tempAAT(m)%AA1=a1(1)(1:1)
        tempAAT(m)%AA3=a2(1)(1:3)
!        PRINT *,m,tempAAT(m)%AA1
      END DO amino
    END IF
  END DO main

  CLOSE (UNIT=oldlognum)

!  PRINT *,'old log closed'

! stupid: if previous run was DNA only, set CFT to E. coli Class II

  IF (.not.oldCodons) THEN
    DO i=1,21
      AAT(i)%NumOfCodons=0
      AAT(i)%NumOfActiveCodons=0
    END DO
    CALL Read_Organism_CFT('ecoli2')
  ELSE

! regenerate CFT and AAT information

    DO i=1,21
      AAT(i)=tempAAT(i)
      AAT(i)%NumOfCodons=0
      AAT(i)%NumOfActiveCodons=0
    END DO

    DO i=1,64
      CFT(i)=tempCFT(i)
      CFT(i)%SeqRC(1:3)=CFT(i)%Seq(1:3)
      CALL RevComplStr(CFT(i)%SeqRC(1:3))
      inner: DO j=1,21
        IF (CFT(i)%AA1.eq.AAT(j)%AA1) THEN
          AAT(j)%NumOfCodons=AAT(j)%NumOfCodons+1
          AAT(j)%Codon(AAT(j)%NumOfCodons)=i
          CFT(i)%AA3=AAT(j)%AA3
          AAT(j)%Freq(AAT(j)%NumOfCodons)=CFT(i)%Freq
!        PRINT *,i,j,AAT(j)%AA1,CFT(i)%seq,CFT(i)%Freq
          EXIT inner
        END IF
      END DO inner
    END DO

! assign numerical values to codon sequence

    DO i=1,64
      DO k=1,3
        SELECT CASE(CFT(i)%Seq(k:k))
          CASE('A')
            CFT(i)%num(k)=-1
          CASE('T')
            CFT(i)%num(k)=1
          CASE('C')
            CFT(i)%num(k)=-3
          CASE('G')
            CFT(i)%num(k)=3
        END SELECT
      END DO
    END DO
  END IF

  CALL Fix_Codons

END SUBROUTINE Read_Old_Codons
SUBROUTINE Read_Old_DNA(str)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: ierr,begin,tempDNAlen,i,j,ne,ns,nlines
  CHARACTER(LEN=100) :: line
  CHARACTER(LEN=22) :: str
  CHARACTER(LEN=9999) :: tempDNAseq
  CHARACTER(LEN=1) :: nt

  IF (TEST0) PRINT *,'Read_Old_DNA'  

  OPEN (UNIT=oldlognum,FILE=oldlogfile,STATUS='old',IOSTAT=ierr)

  main: DO
    READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
    IF (ierr.lt.0) EXIT main   ! dump out if end of file reached
    IF (INDEX(line,str(1:22)).eq.0) CYCLE main ! don't start until DNAstring is found
    begin=1
    tempDNAlen=0
    count: DO
      READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
      IF (ierr.lt.0) EXIT main
      IF (INDEX(line,"The oligonucleotide assembly is:").gt.0) EXIT count
      IF (INDEX(line,"----------").eq.0) THEN
! Just read the line and extract only the DNA characters
        dnaread: DO j=1,LEN_TRIM(line(1:100))
          nt=line(j:j)
          IF (nt.eq.'A'.or.nt.eq.'C'.or.nt.eq.'G'.or.nt.eq.'T'.or.nt.eq.'N'.or. &
              nt.eq.'M'.or.nt.eq.'R'.or.nt.eq.'W'.or.nt.eq.'S'.or.nt.eq.'Y'.or. &
              nt.eq.'K'.or.nt.eq.'B'.or.nt.eq.'D'.or.nt.eq.'H'.or.nt.eq.'V') THEN
            tempDNAlen=tempDNAlen+1
            tempDNAseq(tempDNAlen:tempDNAlen)=nt
          END IF
        END DO dnaread
      END IF
    END DO count
    EXIT main
  END DO main

! make sure the old and new sequences are the same length

  IF (tempDNAlen.ne.DNAlen) THEN

    WRITE(UNIT=console,FMT="(' ')")
    WRITE(UNIT=console,FMT='(16x,"ORIGINAL DNA SEQUENCE")')
    WRITE(UNIT=console,FMT="(1x,a64)") bar64

    nlines=tempDNAlen/60
    DO i=1,nlines
      ne=i*60
      ns=ne-59
      WRITE(UNIT=console,FMT="(i4,1x,a60)") ns,tempDNAseq(ns:ne)
    END DO
    ns=nlines*60+1
    ne=tempDNAlen
    WRITE(UNIT=console,FMT="(i4,1x,a/)") ns,tempDNAseq(ns:ne)

    WRITE(UNIT=console,FMT="(' ')")
    WRITE(UNIT=console,FMT='(16x,"MUTATED DNA SEQUENCE")')
    WRITE(UNIT=console,FMT="(1x,a64)") bar64

    nlines=DNAlen/60
    DO i=1,nlines
      ne=i*60
      ns=ne-59
      WRITE(UNIT=console,FMT="(i4,1x,a60)") ns,CurrDNA%DNAseq(ns:ne)
    END DO
    ns=nlines*60+1
    ne=DNAlen
    WRITE(UNIT=console,FMT="(i4,1x,a/)") ns,CurrDNA%DNAseq(ns:ne)
    CALL Stop_Program("Original and mutated genes are different lengths.")
  END IF

  OLDDNAseq=tempDNAseq

  CLOSE (UNIT=oldlognum)

END SUBROUTINE Read_Old_DNA
LOGICAL FUNCTION Read_Old_Logfile(trial)

! This massive function attempts to read in an old logfile, extract a previous
! trial, and stores the parameters and values in temporary variables.  If
! everything checks out OK, the global variables are then assigned.
!
! It overwrites the following variables:
!
!     CurrDNA%NumOlaps
!     CurrDNA%OlapsPos(,)
!
!     SeqOptimToler
!     OligoConc
!     MgConc
!     SodiumConc
!     OligoLen
!     MeltTemp
!     TBIO
!
!     PTNnum
!     PTN()
!     AAT(21)
!     CFT(64)
!
! It also fills these variables
!
!     OLDDNAseq
!     OLDPROTseq

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER(4) :: ierr,trial,i1,length
  CHARACTER(LEN=25) :: trialString
  CHARACTER(LEN=22) :: DNAString
  CHARACTER(LEN=100) :: line
  CHARACTER(LEN=256) :: a1,a2
  LOGICAL :: found

  found=.FALSE.

  WRITE(line,FMT="(i4)") trial
  line=TRIM(line)
  line=ADJUSTL(line)
  WRITE(trialString,FMT="('PARAMETERS FOR TRIAL ',a)") line(1:4)
  WRITE(DNAstring,FMT="('The DNA sequence #',i4)") trial

  IF (TEST0) PRINT *,"Read_Old_Logfile ",oldlogfile !TEST0

! Do the first assignment of the function

  Read_Old_Logfile=.FALSE.

  IF (trial.gt.0) THEN

    IF (TEST0) PRINT *,"Looking for ",trialString

    OPEN (UNIT=oldlognum,FILE=oldlogfile,STATUS='old',IOSTAT=ierr)
    IF (ierr.ne.0) THEN
      WRITE(UNIT=console,FMT="('')")
      WRITE(UNIT=console,FMT="('Where is the old logfile ""',a,'""?')") oldlogfile(1:LEN_TRIM(oldlogfile))
      CALL Stop_Program("Can't find old logfile")
    ELSE
      main: DO
        READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
        IF (TEST0) PRINT *,"READING : ",line
        
        IF (ierr.lt.0) EXIT main   ! dump out if end of file reached
        length=LEN_TRIM(line)
        i1=INDEX(line,'Job name:')

! Find old jobname and set LogfileOffset

        IF (i1.gt.0) THEN
          READ(line(i1+9:length),*,IOSTAT=ierr) OLDjobname 
          LogfileOffset = i1-1
          IF (TEST0) PRINT *,"LogfileOffset = ",LogfileOffset
          IF (ierr.ne.0) OLDjobname=''
        END IF

        IF (INDEX(line,trialString(1:25)).gt.0) THEN
          found=.TRUE. 
          EXIT main
        END IF
      END DO main
    END IF
    CLOSE (UNIT=oldlognum)

    IF (.not.found) CALL Stop_Program("Can't find chosen trial!")

    CALL Read_Old_DNA(DNAstring)
    CALL Read_Old_Overlaps(DNAString)
    CALL Read_Old_Parameters(trialString)
    CALL Read_Old_Codons
    CALL Read_Old_Patterns
    CALL Read_Old_Protein

    Read_Old_Logfile=.TRUE.

  END IF

END FUNCTION Read_Old_Logfile
SUBROUTINE Read_Old_Overlaps(str)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=22) :: str
  CHARACTER(LEN=100) :: line,line1,line2
  CHARACTER(LEN=1) :: x1,x2
  INTEGER(4) :: i,j,k,n,m,q
  INTEGER(4) :: ierr,linenum,pos,last
  INTEGER(4) :: pt1,pt2,ppt1,ppt2
  INTEGER(4) :: tempNumOlaps,tempOlapsPos(250,2)
  INTEGER,EXTERNAL :: StrToInt
!  REAL,EXTERNAL :: TmCalc

  tempOlapsPos(1,1)=0
  tempNumOlaps=0
  n=1
  m=1
  q=-1
  ppt1=0
  ppt2=0

  IF (TEST0) PRINT *,'Read_Old_Overlaps'  

  OPEN (UNIT=oldlognum,FILE=oldlogfile,STATUS='old',IOSTAT=ierr)

! Find the start position in the file

  main: DO
    READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
    IF (ierr.lt.0) CALL Stop_Program("Can't find overlap block!") 
    IF (INDEX(line,str(1:22)).gt.0) EXIT main 
  END DO main

! calculate expected number of lines

  linenum = (INT((DNAlen-1)/60))+1

! skip ahead

  look: DO
    READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
    IF (ierr.lt.0) CALL Stop_Program("Premature end in overlap block!") 
    IF (INDEX(line,"The oligonucleotide assembly is").gt.0) THEN
      DO i=1,5
        READ(unit=oldlognum,FMT='(A100)',IOSTAT=ierr) line
        IF (ierr.lt.0) CALL Stop_Program("Premature end in overlap block!")
      END DO
      EXIT look
    END IF
  END DO look

! scan the assembly to determine the overlap positions and number of overlaps

  linevi: DO i=1,linenum

! get sense line

    READ(unit=oldlognum,FMT='(A100)',IOSTAT=ierr) line1
    IF (ierr.lt.0) CALL Stop_Program("Premature end in overlap block!") 
    skip: DO

! get antisense line

      READ(unit=oldlognum,FMT='(A100)',IOSTAT=ierr) line2
      IF (ierr.lt.0) CALL Stop_Program("Premature end in overlap block!") 
      IF (INDEX(line2,"*").eq.0) EXIT skip
    END DO skip

! parse the sense and antisense lines

    screen: DO j=1,60
      pos=((i-1)*60)+j

! x1 is the sense strand character at current pos

      x1=line1((j+(LogfileOffset+4)):(j+(LogfileOffset+4)))
      IF ((x1.eq.'A').or.(x1.eq.'C').or.(x1.eq.'G').or.(x1.eq.'T'))  pt1=2
      IF ((x1.eq.'a').or.(x1.eq.'c').or.(x1.eq.'g').or.(x1.eq.'t'))  pt1=1
      IF (x1.eq.' ')                                                 pt1=0

! x2 is the antisense strand character at current pos

      x2=line2((j+(LogfileOffset+4)):(j+(LogfileOffset+4)))
      IF ((x2.eq.'A').or.(x2.eq.'C').or.(x2.eq.'G').or.(x2.eq.'T'))  pt2=2
      IF ((x2.eq.'a').or.(x2.eq.'c').or.(x2.eq.'g').or.(x2.eq.'t'))  pt2=1
      IF (x2.eq.' ')                                                 pt2=0

      IF (tempOlapsPos(1,1).ne.0) THEN
        IF (pt1.ne.ppt1) THEN
          n=((q+1)/2)+1
          IF (ppt1.eq.0) THEN
            tempOlapsPos(m,n) = pos
          ELSE
            tempOlapsPos(m,n) = pos-1
          END IF
          IF (n.eq.2) tempNumOlaps=tempNumOlaps+1
          m=m+((q+1)/2)
          q=q*(-1)
          IF ((pt1.ne.0).and.(ppt1.ne.0)) THEN
            n=((q+1)/2)+1
            tempOlapsPos(m,n) = pos
            IF (n.eq.2) tempNumOlaps=tempNumOlaps+1
            m=m+((q+1)/2)
            q=q*(-1)
          END IF
        END IF
      END IF

      IF (pt2.ne.ppt2) THEN
        n=((q+1)/2)+1
        IF (tempOlapsPos(1,1).eq.0) THEN
          tempOlapsPos(m,n) = pos
        ELSE
          IF (ppt2.eq.0) THEN
            tempOlapsPos(m,n) = pos
          ELSE
            tempOlapsPos(m,n) = pos-1
          END IF
        END IF
        IF (n.eq.2) tempNumOlaps=tempNumOlaps+1
        m=m+((q+1)/2)
        q=q*(-1)
        IF ((pt2.ne.0).and.(ppt2.ne.0)) THEN
          n=((q+1)/2)+1
         tempOlapsPos(m,n) = pos
          IF (n.eq.2) tempNumOlaps=tempNumOlaps+1
          m=m+((q+1)/2)
          q=q*(-1)
        END IF
      END IF

      ppt1=pt1
      ppt2=pt2

    END DO screen

! skip to the next line

    DO j=1,6
      READ(unit=oldlognum,FMT='(A100)',IOSTAT=ierr) line
      IF (ierr.lt.0) CALL Stop_Program("Premature end in overlap block!")
    END DO

  END DO linevi

  CLOSE(UNIT=oldlognum)

  IF (tempNumOlaps.le.0) THEN
    CALL Stop_Program("Zero overlaps in old logfile")
  ELSE
    IF (MOD(tempNumOlaps,2).eq.1) THEN
      DO i=1,tempNumOlaps
        IF (tempOlapsPos(i,1).eq.0.or.tempOlapsPos(i,2).eq.0) THEN
          CALL Stop_Program("Zero old overlaps")
        END IF
        IF (tempOlapsPos(i,1).ge.tempOlapsPos(i,2)) THEN
          CALL Stop_Program("Overlap pos 1 is greater than 2")
        END IF
        IF ((i.gt.1).and.(tempOlapsPos(i,1).le.tempOlapsPos(i-1,2))) THEN
          CALL Stop_Program("Overlap inversion")
        END IF
      END DO
    ELSE
      CALL Stop_Program("Odd number of overlaps in old logfile")
    END IF
  END IF

  CurrDNA%NumOlaps=tempNumOlaps
  DO i=1,tempNumOlaps
    CurrDNA%OlapsPos(i,1)=tempOlapsPos(i,1)
    CurrDNA%OlapsPos(i,2)=tempOlapsPos(i,2)
!    CurrDNA%MeltT(i)=TmCalc(CurrDNA%OlapsPos(i,1),CurrDNA%OlapsPos(i,2))
  END DO

END SUBROUTINE Read_Old_Overlaps
SUBROUTINE Read_Old_Parameters(str)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=100) :: line
  CHARACTER(LEN=25) :: str
  REAL :: tempOligoConc,tempMgConc,tempSodiumConc
  INTEGER(4) :: ierr,tempOligoLen,tempMeltTemp,tempSeqOptimToler
  LOGICAL :: tempTBIO
  REAL,EXTERNAL :: StrToReal
  INTEGER,EXTERNAL :: StrToInt

  IF (TEST0) PRINT *,'Read_Old_Parameters'  

  OPEN (UNIT=oldlognum,FILE=oldlogfile,STATUS='old',IOSTAT=ierr)

  tempTBIO=.FALSE.
  main: DO
    READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
    IF (ierr.lt.0) EXIT main   ! dump out if end of file reached
    IF (INDEX(line,str(1:25)).eq.0) CYCLE main
    parameters: DO
      READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
      IF (ierr.lt.0) EXIT main
      IF (INDEX(line,"DNA sequence").gt.0) EXIT parameters
      SELECT CASE(line((LogfileOffset+15):(LogfileOffset+21)))
        CASE('Oligo S')
          tempOligoLen=StrToInt(line((LogfileOffset+43):(LogfileOffset+46)))
        CASE('Anneali')
          tempMeltTemp=StrToInt(line((LogfileOffset+43):(LogfileOffset+46)))
        CASE('Oligo C')
          tempOligoConc=StrToReal(line((LogfileOffset+43):(LogfileOffset+59)))
        CASE('Sodium ')
          tempSodiumConc=StrToReal(line((LogfileOffset+43):(LogfileOffset+59)))
        CASE('Mg2+ Co')
          tempMgConc=StrToReal(line((LogfileOffset+43):(LogfileOffset+59)))
        CASE('Codon F')
          tempSeqOptimToler=StrToInt(line((LogfileOffset+43):(LogfileOffset+59)))
        CASE('odynami')
          tempTBIO=.TRUE.
        CASE('of even')
          tempOligoLen=StrToInt(line((LogfileOffset+61):(LogfileOffset+64)))
        CASE('-------')
          CYCLE parameters
        CASE('       ')
          CYCLE parameters
      END SELECT
    END DO parameters
  END DO main

  CLOSE (UNIT=oldlognum)

! Check values for robustness

  IF (tempMeltTemp.LT.55.OR.tempMeltTemp.GT.75) THEN
    PRINT *,"MeltTemp bad"
  END IF
  IF (tempOligoLen.LT.20) THEN
    PRINT *,"OligoLen bad"
  END IF

! substitute parameters with old parameters

  OligoConc=tempOligoConc
  MgConc=tempMgConc
  SodiumConc=tempSodiumConc
  OligoLen=tempOligoLen
  MeltTemp=tempMeltTemp
  TBIO=tempTBIO
  SeqOptimToler=tempSeqOptimToler
  MeltTempLo=MeltTemp
  MeltTempHi=MeltTemp
  OligoLenLo=OligoLen
  OligoLenHi=OligoLen

END SUBROUTINE Read_Old_Parameters
SUBROUTINE Read_Old_Patterns

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: ierr,i,j
  CHARACTER(LEN=100) :: line
  TYPE(Pattern) :: RS(250)
  INTEGER(4) :: RSN

  IF (TEST0) PRINT *,'Read_Old_Patterns'  

  OPEN (UNIT=oldlognum,FILE=oldlogfile,STATUS='old',IOSTAT=ierr)

! -- SEQUENCE PATTERNS / RESTRICTION SITES --

  main: DO
    READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
    IF (INDEX(line,"RESTRICTION SITES").gt.0.or.INDEX(line,"SEQUENCE PATTERNS").gt.0) THEN
      READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
!      PRINT *,'restriction sites'
      RSN=0
      restrict: DO
        READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
        IF (INDEX(line,"PARAMETERS FOR TRIAL").gt.0.or.ierr.lt.0) EXIT main
        IF (line((LogfileOffset+1):(LogfileOffset+3)).eq."---".or.line((LogfileOffset+1):(LogfileOffset+3)).eq."   ") CYCLE restrict
        RSN=RSN+1
        RS(RSN)%Name(1:12)=line((LogfileOffset+2):(LogfileOffset+13))
        RS(RSN)%Seq(1:66)=line((LogfileOffset+14):(LogfileOffset+79))
        RS(RSN)%Seq=ADJUSTL(RS(RSN)%Seq)
        RS(RSN)%Seq=TRIM(RS(RSN)%Seq)
        RS(RSN)%Len=LEN_TRIM(RS(RSN)%Seq)
        RS(RSN)%SeqRC=RS(RSN)%Seq
        CALL RevComplStr(RS(RSN)%SeqRC)
        RS(RSN)%SeqRC=ADJUSTL(RS(RSN)%SeqRC)
        RS(RSN)%SeqRC=TRIM(RS(RSN)%SeqRC)

! determine if site is degenerate

        degcheck: DO j=1,RS(RSN)%Len
          IF ((RS(RSN)%Seq(j:j).ne.'A').and.(RS(RSN)%Seq(j:j).ne.'C').and.&
            (RS(RSN)%Seq(j:j).ne.'G').and.(RS(RSN)%Seq(j:j).ne.'T')) THEN
            RS(RSN)%Degen=.TRUE.
            EXIT degcheck
          END IF
        END DO degcheck

! determine if site is self-complementary

        IF (RS(RSN)%Seq(1:RS(RSN)%Len).eq.RS(RSN)%SeqRC(1:RS(RSN)%Len)) THEN
          RS(RSN)%SelfCompl=.TRUE.
        END IF

! Make sure that only unique sites are present

        IF (RSN.gt.1) THEN
          isoschiz: DO j=1,RSN-1
            IF (RS(RSN)%Seq.eq.RS(j)%Seq) RS(RSN)%Isoschiz=.TRUE.
          END DO isoschiz
        END IF

      END DO restrict
    END IF
  END DO main

  DO i=1,RSN
    PTN(i)=RS(i)
  END DO
  PTNnum=RSN

  CLOSE (UNIT=oldlognum)

END SUBROUTINE Read_Old_Patterns
SUBROUTINE Read_Old_Protein

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: ierr,i,j,k,begin,tempPROTlen,nlines,ne,ns
  CHARACTER(LEN=3333) :: tempPROTseq
  CHARACTER(LEN=100) :: line
  INTEGER,EXTERNAL :: StrToInt
  CHARACTER(LEN=1) :: nt
  INTEGER,EXTERNAL :: NT2Int

  IF (TEST0) PRINT *,'Read_Old_Protein'  

  OPEN (UNIT=oldlognum,FILE=oldlogfile,STATUS='old',IOSTAT=ierr)

  tempPROTlen=0
  tempPROTseq=''
  begin=1

  i=0

  main: DO
    READ (UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line
    IF (ierr.lt.0) EXIT main   ! dump out if end of file reached

    IF ((INDEX(line,' PROTEIN ').gt.0).and.(INDEX(line,' SEQUENCE ').gt.0)) THEN
      READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line  ! read bar line
      protein: DO
        READ(UNIT=oldlognum,FMT='(A100)',IOSTAT=ierr) line  ! read sequence
        IF (ierr.lt.0) EXIT main
        IF (INDEX(line,'---------').gt.0) EXIT protein
        IF (INDEX(line,'CODON FREQUENCY').gt.0) EXIT main
        IF (StrToInt(line((LogfileOffset+0):(LogfileOffset+3))).gt.0) THEN
          tempPROTlen = tempPROTlen+LEN_TRIM(line((LogfileOffset+5):(LogfileOffset+64)))
          tempPROTseq(begin:tempPROTlen)=line((LogfileOffset+5):(tempPROTlen-begin+(LogfileOffset+5)))
          begin=tempPROTlen+1
        END IF
      END DO protein
    END IF
  END DO main

! make sure current and old protein sequences are the same length

  IF (tempPROTlen.ne.PROTlen.and.PROTlen.ne.0) THEN

    WRITE(UNIT=console,FMT="(' ')")
    WRITE(UNIT=console,FMT='(16x,"ORIGINAL PROTEIN SEQUENCE")')
    WRITE(UNIT=console,FMT="(1x,a64)") bar64
 
    nlines=tempPROTlen/60
    DO i=1,nlines
      ne=i*60
      ns=ne-59
      WRITE(UNIT=console,FMT="(i4,1x,a60)") ns,tempPROTseq(ns:ne)
    END DO
    ns=nlines*60+1
    ne=tempPROTlen
    WRITE(UNIT=console,FMT="(i4,1x,a/)") ns,tempPROTseq(ns:ne)
 
    WRITE(UNIT=console,FMT="(' ')")
    WRITE(UNIT=console,FMT='(16x,"MUTATED PROTEIN SEQUENCE")')
    WRITE(UNIT=console,FMT="(1x,a64)") bar64
 
    nlines=PROTlen/60
    DO i=1,nlines
      ne=i*60
      ns=ne-59
      WRITE(UNIT=console,FMT="(i4,1x,a60)") ns,PROTseq(ns:ne)
    END DO
    ns=nlines*60+1
    ne=PROTlen
    WRITE(UNIT=console,FMT="(i4,1x,a/)") ns,PROTseq(ns:ne)
 
    CALL Stop_Program("Original and mutated proteins are different lengths.")
  END IF

  IF (tempPROTlen.gt.0) THEN
    OLDPROTseq=tempPROTseq

    mutPROTnum=0

    DO i=1,PROTlen  

! If the protein residues are identical, preserve the old DNA sequence

      IF (PROTseq(i:i).eq.OLDPROTseq(i:i)) THEN
        j=prot2nt(i)
        CurrDNA%DNAseq(j-1:j+1) = OLDDNAseq(j-1:j+1) 

!        PRINT *,'Preserving residue',i,j,OLDDNAseq(j-1:j+1)

! Update prot2cod, nt2cod, and NUMseq arrays

        DO k=1,64
          IF (CurrDNA%DNAseq(j-1:j+1).eq.CFT(k)%Seq) THEN
            CurrDNA%prot2cod(i)=k
            CurrDNA%nt2cod(j)=k
!            PRINT *,'i,j,k',i,j,k,CFT(k)%Seq
          END IF
        END DO
    
        DO k=j-1,j+1
          nt=CurrDNA%DNAseq(k:k)
          IF (nt.eq.'A'.or.nt.eq.'C'.or.nt.eq.'G'.or.nt.eq.'T') CurrDNA%NUMseq(k)=NT2Int(nt)
        END DO

! Otherwise, flag the protein residue for mutation

      ELSE
        mutPROTnum = mutPROTnum+1
        mutPROT2prot(mutPROTnum) = i 
      END IF
    END DO
  END IF

  CLOSE (UNIT=oldlognum)

END SUBROUTINE Read_Old_Protein
SUBROUTINE Read_Organism_CFT(choice)

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=6) :: choice
  INTEGER :: i
  CHARACTER(LEN=12) :: filler="         0.0"

  IF (TEST0) PRINT *,"Read_Organism_CFT" !TEST0

  SELECT CASE(choice(1:6))
    CASE('ecoli2')
      Organism="E. coli Class II"
      DO i=1,64
        CALL Read_Codon_Line(ecoli2CFT(1,i),ecoli2CFT(2,i),filler,ecoli2CFT(3,i))
      END DO
    CASE('ecoli ')
      Organism="E. coli"
      DO i=1,64
        CALL Read_Codon_Line(ecoCFT(1,i),ecoCFT(2,i),filler,ecoCFT(3,i))
      END DO
    CASE('celega')
      Organism="C. elegans"
      DO i=1,64
        CALL Read_Codon_Line(celCFT(1,i),celCFT(2,i),filler,celCFT(3,i))
      END DO
    CASE('drosop')
      Organism="D. melanogaster"
      DO i=1,64
        CALL Read_Codon_Line(dmeCFT(1,i),dmeCFT(2,i),filler,dmeCFT(3,i))
      END DO
    CASE('homosa')
      Organism="H. sapiens"
      DO i=1,64
        CALL Read_Codon_Line(hsaCFT(1,i),hsaCFT(2,i),filler,hsaCFT(3,i))
      END DO
    CASE('musmus')
      Organism="M. musculus"
      DO i=1,64
        CALL Read_Codon_Line(mmuCFT(1,i),mmuCFT(2,i),filler,mmuCFT(3,i))
      END DO
    CASE('pastor')
      Organism="P. pastoris"
      DO i=1,64
        CALL Read_Codon_Line(ppaCFT(1,i),ppaCFT(2,i),filler,ppaCFT(3,i))
      END DO
    CASE('rattus')
      Organism="R. norvegicus"
      DO i=1,64
        CALL Read_Codon_Line(rnoCFT(1,i),rnoCFT(2,i),filler,rnoCFT(3,i))
      END DO
    CASE('saccho')
      Organism="S. cerevesiae"
      DO i=1,64
        CALL Read_Codon_Line(sceCFT(1,i),sceCFT(2,i),filler,sceCFT(3,i))
      END DO
    CASE('xenopu')
      Organism="X. laevis"
      DO i=1,64
        CALL Read_Codon_Line(xlaCFT(1,i),xlaCFT(2,i),filler,xlaCFT(3,i))
      END DO
  END SELECT

END SUBROUTINE Read_Organism_CFT
SUBROUTINE Read_Parameters

! This reads in the parameters from the input file DNAWORKS.inp.  The 
! parameters are parsed using the following keywords and modifiers:

! TBIO
! NOGAPS
! LOGFILE "string" (quotes recommended, but not required)
! TITLE "string" (quotes recommended, but not required)
! EMAIl "string" (quotes recommended, but not required)
! SOLUTIONS #
! MELTING # [LOW #] [HIGh #] [TOLerance #]
! LENGTH # [LOW #] [HIGh #]
! FREQUENCY # [SCORE] [STRICt] [RANDOm]
! REPEAT #
! CONCENTRATION [OLIgo #] [SODium #] [MAGnesium #]
! MISPRIME # [TIP #] [MAX #]
! WEIGHT [AWT #] [CWT #] [GWT #] [HWT #] [LWT #] [MWT #] [PWT #] [RWT #]
! TIMELIMIT #

! The keywords are case insensitive -- they are read from an uppercase only
! line.

! Note: If a string is flanked by quotes (single or double) and is read 
! by free format, the quotes are discarded.  However, to avoid any 
! problems all quoted strings are eliminated when the lines from the input
! file are converted to uppercase.

! Note: All directives must start in the first column against the left margin.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=256) :: a1,a2,a3,a4,a5,a6 ! throw away strings

  INTEGER,EXTERNAL :: StrToInt
  REAL,EXTERNAL :: StrToReal

  INTEGER :: ierr,length,i,j,k,n
  INTEGER :: i1,i2,i3,i4,i5,i6  ! throw away integers

  REAL :: r1,r2,r3,r4,r5,r6  ! throw away reals

  LOGICAL :: skip

  skip = .FALSE.

  IF (TEST0) PRINT *,"Read_Parameters" !TEST0

  main: DO i=1,InputArrayNum

! skip the nucleotide and protein blocks

    IF (INDEX(InputArrayUC(i),'PROT').gt.0) skip=.TRUE.
    IF (INDEX(InputArrayUC(i),'NUCLEOTIDE').eq.1) skip=.TRUE.
    IF (INDEX(InputArrayUC(i),'//').gt.0) skip=.FALSE.
    IF (skip) CYCLE main

! read an individual line, trim it, and change to upper case

    length=LEN_TRIM(InputArray(i))

! Internal testing

    IF (INDEX(InputArrayUC(i),'DNAWORKS INTERNAL TESTING: ').gt.0) THEN
      IF (INDEX(InputArrayUC(i),'TEST0').gt.0) TEST0=.TRUE.
      IF (INDEX(InputArrayUC(i),'TEST1').gt.0) TEST1=.TRUE.
      IF (INDEX(InputArrayUC(i),'TEST2').gt.0) TEST2=.TRUE.
      IF (INDEX(InputArrayUC(i),'TEST3').gt.0) TEST3=.TRUE.
    END IF ! internal testing block

! keyword: TBIO

    IF (INDEX(InputArrayUC(i),'TBIO').eq.1) TBIO=.TRUE. ! tbio block

! keyword: MASCSS

    IF (INDEX(InputArrayUC(i),'MASCSS').gt.0) JACEK=.TRUE. ! jacek block

! keyword: NOGAPS

    IF (INDEX(InputArrayUC(i),'NOGAPS').eq.1) NOGAPS=.TRUE. 

! keyword: QUIET

    IF (INDEX(InputArrayUC(i),'QUIET').gt.0) QUIET=.TRUE. 

! keyword: LOGFILE
!   "<str>" = outputfile

    IF (INDEX(InputArrayUC(i),'LOGFILE').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'LOGFILE')
      IF (i1.gt.0) THEN
        READ(InputArray(i)(i1:length),*,IOSTAT=ierr) a1,a2
        IF (ierr.eq.0) outputfile=a2
      END IF
    END IF ! logfile block

! keyword: TITLE
!   "<str>" = jobname

    IF (INDEX(InputArrayUC(i),'TITLE').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'TITLE')
      IF (i1.gt.0) THEN
        READ(InputArray(i)(i1:length),*,IOSTAT=ierr) a1,a2
        IF (ierr.eq.0) jobname=a2
      END IF
    END IF ! title block

! keyword: MAILPATH
!   "<str>" = MAILPATH

    IF (INDEX(InputArrayUC(i),'MAILPATH').gt.0) THEN
      i1=INDEX(InputArrayUC(i),'MAILPATH')
      IF (i1.gt.0) THEN
        READ(InputArray(i)(i1:length),*,IOSTAT=ierr) a1,a2
        IF (ierr.eq.0) MAILPATH=a2
      END IF
    END IF ! MAILPATH block

! keyword: EMAIl
!   "<str>" = email

    IF (INDEX(InputArrayUC(i),'EMAI').gt.0) THEN
      i1=INDEX(InputArrayUC(i),'EMAI')
      IF (i1.gt.0) THEN
        READ(InputArray(i)(i1:length),*,IOSTAT=ierr) a1,a2
        IF (ierr.eq.0) THEN
          i3=INDEX(a2,"@")
          i4=INDEX(a2,"@",.TRUE.)
          IF (i3.ne.0.and.i3.eq.i4) email=a2
        END IF
      END IF
    END IF ! email block

! keyword: SOLUTIONS
!   # = NumberOfSolutions

    IF (INDEX(InputArrayUC(i),'SOLUTIONS').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'SOLUTIONS')
      IF (i1.gt.0)  READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i2
      IF (ierr.eq.0) THEN
        NumberOfSolutions=i2
        IF (NumberOfSolutions.gt.99) NumberOfSolutions=99
        IF (NumberOfSolutions.lt.1) NumberOfSolutions=1
      END IF
    END IF

! keyword: TIMELIMIT
!   # = MainTimeLimit

    IF (INDEX(InputArrayUC(i),'TIMELIMIT').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'TIMELIMIT')
      IF (i1.gt.0)  READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i2
      IF (ierr.eq.0) THEN
        MainTimeLimit=i2
        IF (MainTimeLimit.lt.1) MainTimeLimit=0
      END IF
    END IF

! keyword: MELTING
!   # = MeltTemp,MeltTempHi,MeltTempLo
!   HIG # = MeltTempHi      Ideal melting temperature (upper limit)
!   LOW # = MeltTempLo      Ideal melting temperature (lower limit)
!   TOL # = MeltTol         Tolerance for melting temperature deviation

    IF (INDEX(InputArrayUC(i),'MELTING').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'MELTING')
      IF (i1.gt.0) THEN
        READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i2
        IF (ierr.eq.0) THEN                             ! single number
          MeltTemp=i2
          MeltTempLo=i2
          MeltTempHi=i2
        ELSE
          i1=INDEX(InputArrayUC(i),'LOW')
          i2=INDEX(InputArrayUC(i),'HIG')
          i3=INDEX(InputArrayUC(i),'TOL')
          IF (i1.gt.0) THEN                             ! LOW defined
            READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i4
            IF (ierr.eq.0) MeltTempLo=i4
          END IF
          IF (i2.gt.0) THEN                             ! HIGh defined
            READ(InputArrayUC(i)(i2:length),*,IOSTAT=ierr) a1,i4
            IF (ierr.eq.0) MeltTempHi=i4
          END IF
          IF (i3.gt.0) THEN                             ! HIGh defined
            READ(InputArrayUC(i)(i3:length),*,IOSTAT=ierr) a1,i4
            IF (ierr.eq.0) MeltTol=i4
          END IF
        END IF

! LOW defined, HIGh not defined

        IF (i1.gt.0.and.i2.eq.0) MeltTempHi = MeltTempLo

! HIGh defined, LOW not defined

        IF (i2.gt.0.and.i1.eq.0) MeltTempLo = MeltTempHi

! Imbalanced Hi and Lo MeltTemps

        IF (MeltTempHi.lt.MeltTempLo) CALL IntSwap(MeltTempHi,MeltTempLo)

      END IF
    END IF ! melt block

! keyword: LENGTH
!   # = OligoLen,OligoLenHi,OligoLenLo
!   HIG # = OligoLenHi      user input oligo size (upper limit)
!   LOW # = OligoLenLo      user input oligo size (lower limit)
!   RANDO = OligoLenRandom  allow oligolen to vary between 20 and chosen length

    IF (INDEX(InputArrayUC(i),'LENGTH').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'LENGTH')
      IF (i1.gt.0) THEN
        IF (INDEX(InputArrayUC(i),'RANDO').gt.0) OligoLenRandom=.TRUE.
        READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i2
        IF (ierr.eq.0) THEN                             ! single number
          OligoLen=i2 
          OligoLenLo=i2
          OligoLenHi=i2
        ELSE
          i1=INDEX(InputArrayUC(i),'LOW')
          i2=INDEX(InputArrayUC(i),'HIG')
          IF (i1.gt.0) THEN                             ! LOW defined
            READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i4
            IF (ierr.eq.0) OligoLenLo=i4
          END IF
          IF (i2.gt.0) THEN                             ! HIGh defined
            READ(InputArrayUC(i)(i2:length),*,IOSTAT=ierr) a1,i4
            IF (ierr.eq.0) OligoLenHi=i4
          END IF
        END IF

! LOW defined, HIGh not defined

        IF (i1.gt.0.and.i2.eq.0) OligoLenHi = OligoLenLo

! HIGh defined, LOW not defined

        IF (i2.gt.0.and.i1.eq.0) OligoLenLo = OligoLenHi

! Imbalanced Hi and Lo MeltTemps
        IF (OligoLenHi.lt.OligoLenLo) CALL IntSwap(OligoLenHi,OligoLenLo)
      END IF
    END IF ! length block

! keyword: FREQUENCY 
!   THR # = SeqOptimToler
!   SCORE = ScoreCodons
!   STRIC = CodonStrict
!   RANDO = CodonRandom

    IF (INDEX(InputArrayUC(i),'FREQUENCY').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'FREQUENCY')
      IF (i1.gt.0) THEN
        IF (INDEX(InputArrayUC(i),'STRIC').gt.0) CodonStrict=.TRUE.
        IF (INDEX(InputArrayUC(i),'SCORE').gt.0) ScoreCodons=.TRUE.
        IF (INDEX(InputArrayUC(i),'RANDO').gt.0) CodonRandom=.TRUE.
        i1=INDEX(InputArrayUC(i),'THR')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i2
          IF (ierr.eq.0) THEN
            SeqOptimToler=i2
            IF (SeqOptimToler.gt.100) SeqOptimToler=100
            IF (SeqOptimToler.lt.0) SeqOptimToler=0
          END IF
        END IF
      END IF
    END IF ! frequency block

! keyword: REPEAT
!   # = RepLen

    IF (INDEX(InputArrayUC(i),'REPEAT').gt.0) THEN
      i1=INDEX(InputArrayUC(i),'REPEAT')
      IF (i1.gt.0) THEN
        READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i2
        IF (ierr.eq.0) RepLen=i2
      END IF
    END IF ! repeat block

! keyword: CONCENTRATION
!   OLI # = OligoConcTemp   in nM
!   SOD # = SodiumConc      in mM
!   MAG # = MgConc          in mM

    IF (INDEX(InputArrayUC(i),'CONCENTRATION').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'CONCENTRATION')
      IF (i1.gt.0) THEN
        i1=INDEX(InputArrayUC(i),'OLI')
        i2=INDEX(InputArrayUC(i),'SOD')
        i3=INDEX(InputArrayUC(i),'MAG')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) OligoConc=r1
        END IF
        IF (i2.gt.0) THEN
          READ(InputArrayUC(i)(i2:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) SodiumConc=r1
        END IF
        IF (i3.gt.0) THEN
          READ(InputArrayUC(i)(i3:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) MgConc=r1
        END IF
      END IF
    END IF ! concentration block

! keyword: MISPRIME
!     # = MPLn            length of  misprimes
!   MAX = MaxNonId        maximum number of non-identical nts in misprime
!   TIP = MPTip           identical tip of the misprime, in nts

    IF (INDEX(InputArrayUC(i),'MISPRIME').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'MISPRIME')
      IF (i1.gt.0) THEN
        READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,i1
        IF (i1.gt.0) MPLn=i1
  
        i2=INDEX(InputArrayUC(i),'TIP')
        i3=INDEX(InputArrayUC(i),'MAX')
  
        IF (i2.gt.0) THEN
          READ(InputArrayUC(i)(i2:length),*,IOSTAT=ierr) a1,i2
          IF (ierr.eq.0) MPTip=i2
        END IF
        IF (i3.gt.0) THEN
          READ(InputArrayUC(i)(i3:length),*,IOSTAT=ierr) a1,i3
          IF (ierr.eq.0) MaxNonId=i3
        END IF
      END IF
    END IF ! misprime block

! keyword: WEIGHT
!   TWT # = Twt ! weight for MeltTm scoring
!   CWT # = Cwt ! weight for codon scoring
!   RWT # = Rwt ! weight for repeat scoring
!   MWT # = Mwt ! weight for mispriming scoring
!   GWT # = Gwt ! weight for GC scoring
!   AWT # = Awt ! weight for AT scoring
!   LWT # = Lwt ! weight for length scoring
!   PWT # = Pwt ! weight for pattern scoring
!   FWT # = Fwt ! weight for gap fixing

    IF (INDEX(InputArrayUC(i),'WEIGHT').eq.1) THEN
      i1=INDEX(InputArrayUC(i),'WEIGHT')
      IF (i1.gt.0) THEN
        i1=INDEX(InputArrayUC(i),'FWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Fwt=r1
        END IF
  
        i1=INDEX(InputArrayUC(i),'TWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Twt=r1
        END IF
    
        i1=INDEX(InputArrayUC(i),'CWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Cwt=r1
        END IF
    
        i1=INDEX(InputArrayUC(i),'RWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Rwt=r1
        END IF
    
        i1=INDEX(InputArrayUC(i),'MWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Mwt=r1
        END IF
    
        i1=INDEX(InputArrayUC(i),'GWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Gwt=r1
        END IF
  
        i1=INDEX(InputArrayUC(i),'AWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Awt=r1
        END IF
    
        i1=INDEX(InputArrayUC(i),'LWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Lwt=r1
        END IF
  
        i1=INDEX(InputArrayUC(i),'PWT')
        IF (i1.gt.0) THEN
          READ(InputArrayUC(i)(i1:length),*,IOSTAT=ierr) a1,r1
          IF (ierr.eq.0) Pwt=r1
        END IF
      END IF
    END IF ! weight block

  END DO main

END SUBROUTINE Read_Parameters
SUBROUTINE Read_Pattern_Line(str,seq)
!
! This subroutine assigns the pattern type

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j
  CHARACTER(LEN=80) :: str,seq

  IF (TEST1) PRINT *,"Read_Pattern_Line" !TEST0

  CALL ToUpperCase(seq)

  IF (str.ne.''.and.seq.ne.'') THEN
    PTNnum=PTNnum+1
    i=PTNnum
  
    PTN(i)%Name=ADJUSTL(str)
    PTN(i)%Seq=ADJUSTL(seq)
    PTN(i)%Seq=TRIM(PTN(i)%Seq)
    PTN(i)%Len=LEN_TRIM(PTN(i)%Seq)
  
    PTN(i)%SeqRC=PTN(i)%Seq
    CALL RevComplStr(PTN(i)%SeqRC)
  
    PTN(i)%SeqRC=ADJUSTL(PTN(i)%SeqRC)
    PTN(i)%SeqRC=TRIM(PTN(i)%SeqRC)
  
    degcheck: DO j=1,PTN(i)%Len
      IF ((PTN(i)%Seq(j:j).ne.'A').and.(PTN(i)%Seq(j:j).ne.'C').and.&
          (PTN(i)%Seq(j:j).ne.'G').and.(PTN(i)%Seq(j:j).ne.'T')) THEN
        PTN(i)%Degen=.TRUE.
        EXIT degcheck
      END IF
    END DO degcheck
  
    IF (PTN(i)%Seq(1:PTN(i)%Len).eq.PTN(i)%SeqRC(1:PTN(i)%Len)) THEN
      PTN(i)%SelfCompl=.TRUE.
    END IF

! Make sure that only unique sites are present

    IF (PTNnum.gt.1) THEN
      isoschiz: DO j=1,PTNnum-1
        IF (PTN(i)%Seq.eq.PTN(j)%Seq) PTN(i)%Isoschiz=.TRUE.
      END DO isoschiz
    END IF
  END IF

END SUBROUTINE Read_Pattern_Line
SUBROUTINE Read_Patterns

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=256) :: a1,a2,a3,a4,a5,a6 ! throw away strings

  INTEGER,EXTERNAL :: StrToInt
  REAL,EXTERNAL :: StrToReal

  INTEGER :: ierr,i,j,i1,i2,start,finish

  IF (TEST0) PRINT *,"Read_Patterns" !TEST0

  main: DO i=1,InputArrayNum

! read an individual line, trim it, and change to upper case

! PATTERN

    IF (INDEX(InputArrayUC(i),'PATTERN').eq.1) THEN

! Find where pattern block ends

      start=i+1
      seek: DO j=i,InputArrayNum
        IF (INDEX(InputArrayUC(j),'//').ne.0) THEN
          finish=j
          EXIT seek
        END IF
      END DO seek
    
      pisser: DO j=start,finish
        a1=""
        a2=""
        READ(InputArray(j),*,IOSTAT=ierr) a1,a2
        IF (ierr.ne.0) CYCLE pisser
        CALL Read_Pattern_Line(a1,a2)
      END DO pisser
    END IF

  END DO main

END SUBROUTINE Read_Patterns
