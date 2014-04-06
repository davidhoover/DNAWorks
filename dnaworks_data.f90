MODULE dnaworks_data

  IMPLICIT NONE
  SAVE

! GLOBAL

  INTEGER :: console=6            ! print to console
  INTEGER :: inputnum=9           ! input files
  INTEGER :: outputnum=10         ! output logfile
  INTEGER :: oldlognum=11         ! old logfile output 
  INTEGER :: PrevTrial=0          ! previous trial to fix oligos
  INTEGER :: OligoLen=40          ! user input oligo size
  INTEGER :: OligoLenHi=40        ! user input oligo size (upper limit)
  INTEGER :: OligoLenLo=40        ! user input oligo size (lower limit)
  LOGICAL :: OligoLenRandom=.FALSE.     ! allow oligolen to vary between 20
  INTEGER :: MeltTemp=60          ! Ideal melting temperature
  INTEGER :: MeltTempHi=60        ! Ideal melting temperature (upper limit)
  INTEGER :: MeltTempLo=60        ! Ideal melting temperature (lower limit)
  INTEGER :: MeltTol=1            ! Tolerance for melting temperature deviation
  INTEGER :: SeqOptimToler=50     ! Lowest allowed codon frequency
  INTEGER :: TotalNumberOfSolutions
  INTEGER :: NumberOfSolutions=1
  INTEGER :: RepLen=8             ! determines the size of repeats to minimize
  INTEGER :: MPLn=18              ! length of  misprimes
  INTEGER :: MPTip=6              ! identical tip of the misprime, in nts
  INTEGER :: MaxPROTlen=3333      ! maximum number of protein residues
  INTEGER :: MaxDNAlen=9999       ! maximum number of nucleotide residues
  INTEGER :: MaxNonId=8           ! maximum number of non-identical nts in misprime
  INTEGER :: MutProtPos=0         ! which codon should be mutated
  INTEGER :: MutNtPos(3)          ! which nts are mutated
  INTEGER :: MutNtNum=0           ! how many nts are mutated (zero if none)
  INTEGER :: nt2aa(9999)          ! DNApos to aa (1-21) or 0
  INTEGER :: nt2overlap(9999)     ! DNApos to overlap or 0
  INTEGER :: nt2Solig(9999)       ! DNApos to oligo or 0
  INTEGER :: nt2Aolig(9999)       ! DNApos to antisense oligo or 0
  INTEGER :: nt2prot(9999)        ! DNApos to PROTpos or 0
  INTEGER :: prot2aa(3333)        ! PROTpos to aa (1-21)
  INTEGER :: prot2nt(3333)        ! PROTpos to DNApos (middle nt of codon)
  INTEGER :: DNAlen=0             ! the length of the entire DNA
  INTEGER :: PROTlen=0            ! the length of the all the proteins
  INTEGER :: mutPROTnum=0         ! the number of mutated aa
  INTEGER :: mutPROT2prot(3333)   ! mutPROTpos to PROTpos
  INTEGER :: NumberOfChains=0     ! number of isolated protein chains
  INTEGER :: prot2chain(3333)     ! prot pos to protein number (NumberOfChains)
  LOGICAL :: ChainReverse(99)     ! true if chain is reversed, indexed by
  LOGICAL :: ChainGapFix(99)      ! true if chain is reversed, indexed by

  CHARACTER(LEN=9999) :: INITseq=''     ! initial input sequence (DNA and prot)

! for degenerate nt
  CHARACTER(LEN=9999) :: ORIGDNAseq=''  ! the original dna sequence
  INTEGER :: NumDegPos=0           ! total number of degenerate nt
  INTEGER :: DegPos(999)           ! positions of degenerate sequences
  INTEGER :: CurrSolutionNo=0

  INTEGER :: INITlen=0             ! the length of the initial input sequence
  INTEGER :: NumberOfSeq=0         ! number of sequences, DNA or protein
  INTEGER :: INIT2Seq(9999)        ! 
  LOGICAL :: SeqIsProt(99)         ! true if sequence is prot, false if DNA
  LOGICAL :: SeqReverse(99)        ! true if sequence is reversed, indexed by seq number (NumberOfSeq)
  LOGICAL :: SeqGapFix(99)         ! true if sequence is to be gapfixed, indexed by seq number (NumberOfSeq)
  CHARACTER(LEN=80) :: email=''
  CHARACTER(LEN=80) :: jobname=''
  CHARACTER(LEN=80) :: OLDjobname=''
  CHARACTER(LEN=30) :: oldlogfile="OLDLOGFILE.txt"
  CHARACTER(LEN=30) :: inputfile="DNAWORKS.inp"
  CHARACTER(LEN=30) :: outputfile="LOGFILE.txt"
  CHARACTER(LEN=256) :: InputArray(9999) ! contents of DNAWORKS.inp
  CHARACTER(LEN=256) :: InputArrayUC(9999) ! contents of DNAWORKS.inp, uppercase
  INTEGER :: InputArrayNum              ! number of lines in DNAWORKS.inp

  CHARACTER(LEN=9999) :: SCRATCH=''     ! scratch string for various calls
  CHARACTER(LEN=9999) :: OLDDNAseq=''   ! DNA sequence from previous trial
  CHARACTER(LEN=3333) :: PROTseq=''      ! protein sequence
  CHARACTER(LEN=3333) :: OLDPROTseq=''   ! protein sequence from previous trial
  
  CHARACTER(LEN=64) :: bar64 = "----------------------------------------------------------------"
  CHARACTER(LEN=80) :: bar80 = "--------------------------------------------------------------------------------"
  INTEGER :: MainTimeLimit=0       ! time limit for entire run
  INTEGER :: MainTimeStart              ! for time control of the program
  REAL :: Twt=1.0                       ! weight for MeltTm scoring
  REAL :: Cwt=1.0                       ! weight for codon scoring
  REAL :: Rwt=1.0                       ! weight for repeat scoring
  REAL :: Mwt=1.0                       ! weight for mispriming scoring
  REAL :: Gwt=1.0                       ! weight for GC scoring
  REAL :: Awt=1.0                       ! weight for AT scoring
  REAL :: Lwt=1.0                       ! weight for length scoring
  REAL :: Pwt=1.0                       ! weight for pattern scoring
  REAL :: Fwt=1.0                       ! weight for gap fixing
  REAL :: XScore(3333)                  ! cocon-based total score for mutation
  LOGICAL :: CodonStrict=.FALSE.        ! use strict frequency threshold
  LOGICAL :: ScoreCodons=.FALSE.         ! calculate codon scores
  LOGICAL :: CodonRandom=.FALSE.        ! translate using random codons
  LOGICAL :: MutantRun=.FALSE.          ! if this is a mutation only run
  LOGICAL :: GapFix=.FALSE.             ! are any positions fixed in the gaps?
  INTEGER :: LogfileOffset=1            ! how many blank characters precede the line?
  LOGICAL :: JACEK=.FALSE.
  CHARACTER(LEN=80) :: MAILPATH="/usr/bin/Mail"
  LOGICAL :: TBIO=.FALSE.
  LOGICAL :: NOGAPS=.FALSE.             ! if no gaps are desired
  LOGICAL :: QUIET=.FALSE.
  LOGICAL :: FAST=.FALSE.               ! cut corners
  LOGICAL :: TimesUp=.FALSE.
  LOGICAL :: SequenceTranslated=.FALSE.   ! if false, generate all scores; if
                                        ! true, only generate scores that
                                        ! change when overlaps change
  REAL :: OligoConc=2e-7                ! 200 nM oligo
  REAL :: SodiumConc=5e-2               ! 50 mM sodium
  REAL :: MgConc=2e-3                   ! 2 mM magnesium
  REAL :: RGasConstant=1.9872           ! gas constant
  REAL :: Kelvin=273.15                 ! conversion from kelvin to celsius
  REAL :: OligoCorr                     ! correction factor for oligo conc.
  REAL :: OligoCorrSC                   ! correction factor for self-comp. oligo
  REAL :: SaltCorr                      ! correction factor for cations

! PATTERNS

  TYPE Pattern
    CHARACTER(LEN=80) :: SeqRC
    CHARACTER(LEN=80) :: Seq
    INTEGER :: Len
    CHARACTER(LEN=80) :: Name
    LOGICAL :: SelfCompl
    LOGICAL :: Degen
    LOGICAL :: Isoschiz
  END TYPE

  TYPE(Pattern) :: PTN(999)
  INTEGER :: PTNnum=0

! SOLUTIONS

  TYPE Tally
    REAL :: InitScore
    REAL :: FinaScore
    REAL :: TmRange
    INTEGER :: NumOligs
    INTEGER :: Oligo
    INTEGER :: MeltT
    INTEGER :: LongestOligo
    INTEGER :: Repeats
    INTEGER :: Misprimes
    INTEGER :: LowestOlap
  END TYPE

  TYPE(Tally) :: FinalScore(9999)

! TEST

  TYPE Test_Tally
    REAL :: Score
    INTEGER :: Oligo
    INTEGER :: MeltT
    INTEGER :: Count
    INTEGER :: Time
  END TYPE

  TYPE(Test_Tally) :: Test_Scores(400)

! TABLES

  TYPE KnownCodon
    CHARACTER(LEN=3) :: Seq
    CHARACTER(LEN=3) :: AA3
    CHARACTER(LEN=1) :: AA1
    CHARACTER(LEN=3) :: SeqRC    ! Reverse complement of sequence
    INTEGER :: num(3)         ! numerical representation of codon
    INTEGER :: numRC(3)         ! numerical representation of codon
    REAL :: Freq
    REAL :: Number
    LOGICAL :: Check
  END TYPE KnownCodon

  TYPE(KnownCodon) :: CFT(64)  ! Codon Frequency Table

  TYPE KnownAA
    CHARACTER(LEN=3) :: AA3
    CHARACTER(LEN=1) :: AA1
    REAL :: Freq(10)
    REAL :: NumberSum
    INTEGER :: NumOfCodons
    INTEGER :: NumOfActiveCodons
    INTEGER :: Codon(10)
  END TYPE KnownAA

  TYPE(KnownAA) :: AAT(21)     ! Amino Acid Table

! Degenerate sequences

  TYPE DegenerateSeq           ! Table of degenerate sequences
    CHARACTER(LEN=1) :: DegNT
    INTEGER :: NumOfNT
    INTEGER :: NumSeq(4)
    CHARACTER(LEN=1) :: Seq(4)
  END TYPE DegenerateSeq

  TYPE(DegenerateSeq) :: DegenSeq(11)

! PRE-EXISTANT CFTs

  CHARACTER(LEN=30) :: Organism

  CHARACTER(LEN=5),DIMENSION(3,64) :: ecoli2CFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.044","Gly  ","GGA  ","0.020","Gly  ","GGT  ","0.508","Gly  ","GGC  ","0.428",&
"Glu  ","GAG  ","0.247","Glu  ","GAA  ","0.754","Asp  ","GAT  ","0.461","Asp  ","GAC  ","0.540",&
"Val  ","GTG  ","0.268","Val  ","GTA  ","0.200","Val  ","GTT  ","0.398","Val  ","GTC  ","0.135",&
"Ala  ","GCG  ","0.323","Ala  ","GCA  ","0.240","Ala  ","GCT  ","0.275","Ala  ","GCC  ","0.161",&
"Arg  ","AGG  ","0.003","Arg  ","AGA  ","0.006","Ser  ","AGT  ","0.045","Ser  ","AGC  ","0.243",&
"Lys  ","AAG  ","0.215","Lys  ","AAA  ","0.786","Asn  ","AAT  ","0.173","Asn  ","AAC  ","0.828",&
"Met  ","ATG  ","1.000","Ile  ","ATA  ","0.006","Ile  ","ATT  ","0.335","Ile  ","ATC  ","0.659",&
"Thr  ","ACG  ","0.127","Thr  ","ACA  ","0.047","Thr  ","ACT  ","0.291","Thr  ","ACC  ","0.536",&
"Trp  ","TGG  ","1.000","End  ","TGA  ","0.352","Cys  ","TGT  ","0.389","Cys  ","TGC  ","0.612",&
"End  ","TAG  ","0.076","End  ","TAA  ","0.630","Tyr  ","TAT  ","0.352","Tyr  ","TAC  ","0.648",&
"Leu  ","TTG  ","0.055","Leu  ","TTA  ","0.034","Phe  ","TTT  ","0.291","Phe  ","TTC  ","0.709",&
"Ser  ","TCG  ","0.074","Ser  ","TCA  ","0.048","Ser  ","TCT  ","0.324","Ser  ","TCC  ","0.266",&
"Arg  ","CGG  ","0.008","Arg  ","CGA  ","0.011","Arg  ","CGT  ","0.643","Arg  ","CGC  ","0.330",&
"Gln  ","CAG  ","0.814","Gln  ","CAA  ","0.187","His  ","CAT  ","0.298","His  ","CAC  ","0.702",&
"Leu  ","CTG  ","0.767","Leu  ","CTA  ","0.008","Leu  ","CTT  ","0.056","Leu  ","CTC  ","0.080",&
"Pro  ","CCG  ","0.719","Pro  ","CCA  ","0.153","Pro  ","CCT  ","0.112","Pro  ","CCC  ","0.016"/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: celCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.08 ","Gly  ","GGA  ","0.59 ","Gly  ","GGT  ","0.20 ","Gly  ","GGC  ","0.12 ", &
"Glu  ","GAG  ","0.38 ","Glu  ","GAA  ","0.62 ","Asp  ","GAT  ","0.68 ","Asp  ","GAC  ","0.32 ", &
"Val  ","GTG  ","0.23 ","Val  ","GTA  ","0.16 ","Val  ","GTT  ","0.39 ","Val  ","GTC  ","0.22 ", &
"Ala  ","GCG  ","0.13 ","Ala  ","GCA  ","0.31 ","Ala  ","GCT  ","0.36 ","Ala  ","GCC  ","0.20 ", &
"Arg  ","AGG  ","0.08 ","Arg  ","AGA  ","0.29 ","Ser  ","AGT  ","0.15 ","Ser  ","AGC  ","0.10 ", &
"Lys  ","AAG  ","0.41 ","Lys  ","AAA  ","0.59 ","Asn  ","AAT  ","0.62 ","Asn  ","AAC  ","0.38 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.16 ","Ile  ","ATT  ","0.53 ","Ile  ","ATC  ","0.31 ", &
"Thr  ","ACG  ","0.15 ","Thr  ","ACA  ","0.34 ","Thr  ","ACT  ","0.32 ","Thr  ","ACC  ","0.18 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.39 ","Cys  ","TGT  ","0.55 ","Cys  ","TGC  ","0.45 ", &
"End  ","TAG  ","0.18 ","End  ","TAA  ","0.44 ","Tyr  ","TAT  ","0.56 ","Tyr  ","TAC  ","0.44 ", &
"Leu  ","TTG  ","0.23 ","Leu  ","TTA  ","0.11 ","Phe  ","TTT  ","0.49 ","Phe  ","TTC  ","0.51 ", &
"Ser  ","TCG  ","0.15 ","Ser  ","TCA  ","0.25 ","Ser  ","TCT  ","0.21 ","Ser  ","TCC  ","0.13 ", &
"Arg  ","CGG  ","0.09 ","Arg  ","CGA  ","0.23 ","Arg  ","CGT  ","0.21 ","Arg  ","CGC  ","0.10 ", &
"Gln  ","CAG  ","0.34 ","Gln  ","CAA  ","0.66 ","His  ","CAT  ","0.60 ","His  ","CAC  ","0.40 ", &
"Leu  ","CTG  ","0.14 ","Leu  ","CTA  ","0.09 ","Leu  ","CTT  ","0.25 ","Leu  ","CTC  ","0.17 ", &
"Pro  ","CCG  ","0.20 ","Pro  ","CCA  ","0.53 ","Pro  ","CCT  ","0.18 ","Pro  ","CCC  ","0.09 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: dmeCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.07 ","Gly  ","GGA  ","0.28 ","Gly  ","GGT  ","0.21 ","Gly  ","GGC  ","0.43 ", &
"Glu  ","GAG  ","0.67 ","Glu  ","GAA  ","0.33 ","Asp  ","GAT  ","0.53 ","Asp  ","GAC  ","0.47 ", &
"Val  ","GTG  ","0.47 ","Val  ","GTA  ","0.11 ","Val  ","GTT  ","0.18 ","Val  ","GTC  ","0.24 ", &
"Ala  ","GCG  ","0.19 ","Ala  ","GCA  ","0.17 ","Ala  ","GCT  ","0.19 ","Ala  ","GCC  ","0.45 ", &
"Arg  ","AGG  ","0.11 ","Arg  ","AGA  ","0.09 ","Ser  ","AGT  ","0.14 ","Ser  ","AGC  ","0.25 ", &
"Lys  ","AAG  ","0.70 ","Lys  ","AAA  ","0.30 ","Asn  ","AAT  ","0.44 ","Asn  ","AAC  ","0.56 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.19 ","Ile  ","ATT  ","0.34 ","Ile  ","ATC  ","0.47 ", &
"Thr  ","ACG  ","0.26 ","Thr  ","ACA  ","0.20 ","Thr  ","ACT  ","0.17 ","Thr  ","ACC  ","0.38 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.25 ","Cys  ","TGT  ","0.29 ","Cys  ","TGC  ","0.71 ", &
"End  ","TAG  ","0.33 ","End  ","TAA  ","0.41 ","Tyr  ","TAT  ","0.37 ","Tyr  ","TAC  ","0.63 ", &
"Leu  ","TTG  ","0.18 ","Leu  ","TTA  ","0.05 ","Phe  ","TTT  ","0.37 ","Phe  ","TTC  ","0.63 ", &
"Ser  ","TCG  ","0.20 ","Ser  ","TCA  ","0.09 ","Ser  ","TCT  ","0.08 ","Ser  ","TCC  ","0.23 ", &
"Arg  ","CGG  ","0.15 ","Arg  ","CGA  ","0.15 ","Arg  ","CGT  ","0.16 ","Arg  ","CGC  ","0.33 ", &
"Gln  ","CAG  ","0.70 ","Gln  ","CAA  ","0.30 ","His  ","CAT  ","0.40 ","His  ","CAC  ","0.60 ", &
"Leu  ","CTG  ","0.43 ","Leu  ","CTA  ","0.09 ","Leu  ","CTT  ","0.10 ","Leu  ","CTC  ","0.15 ", &
"Pro  ","CCG  ","0.29 ","Pro  ","CCA  ","0.25 ","Pro  ","CCT  ","0.13 ","Pro  ","CCC  ","0.33 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: hsaCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.25 ","Gly  ","GGA  ","0.25 ","Gly  ","GGT  ","0.16 ","Gly  ","GGC  ","0.34 ", &
"Glu  ","GAG  ","0.58 ","Glu  ","GAA  ","0.42 ","Asp  ","GAT  ","0.46 ","Asp  ","GAC  ","0.54 ", &
"Val  ","GTG  ","0.47 ","Val  ","GTA  ","0.12 ","Val  ","GTT  ","0.18 ","Val  ","GTC  ","0.24 ", &
"Ala  ","GCG  ","0.11 ","Ala  ","GCA  ","0.23 ","Ala  ","GCT  ","0.26 ","Ala  ","GCC  ","0.40 ", &
"Arg  ","AGG  ","0.21 ","Arg  ","AGA  ","0.21 ","Ser  ","AGT  ","0.15 ","Ser  ","AGC  ","0.24 ", &
"Lys  ","AAG  ","0.57 ","Lys  ","AAA  ","0.43 ","Asn  ","AAT  ","0.47 ","Asn  ","AAC  ","0.53 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.17 ","Ile  ","ATT  ","0.36 ","Ile  ","ATC  ","0.47 ", &
"Thr  ","ACG  ","0.11 ","Thr  ","ACA  ","0.28 ","Thr  ","ACT  ","0.25 ","Thr  ","ACC  ","0.36 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.47 ","Cys  ","TGT  ","0.45 ","Cys  ","TGC  ","0.55 ", &
"End  ","TAG  ","0.23 ","End  ","TAA  ","0.30 ","Tyr  ","TAT  ","0.44 ","Tyr  ","TAC  ","0.56 ", &
"Leu  ","TTG  ","0.13 ","Leu  ","TTA  ","0.08 ","Phe  ","TTT  ","0.46 ","Phe  ","TTC  ","0.54 ", &
"Ser  ","TCG  ","0.06 ","Ser  ","TCA  ","0.15 ","Ser  ","TCT  ","0.19 ","Ser  ","TCC  ","0.22 ", &
"Arg  ","CGG  ","0.20 ","Arg  ","CGA  ","0.11 ","Arg  ","CGT  ","0.08 ","Arg  ","CGC  ","0.19 ", &
"Gln  ","CAG  ","0.74 ","Gln  ","CAA  ","0.26 ","His  ","CAT  ","0.42 ","His  ","CAC  ","0.58 ", &
"Leu  ","CTG  ","0.40 ","Leu  ","CTA  ","0.07 ","Leu  ","CTT  ","0.13 ","Leu  ","CTC  ","0.20 ", &
"Pro  ","CCG  ","0.11 ","Pro  ","CCA  ","0.28 ","Pro  ","CCT  ","0.28 ","Pro  ","CCC  ","0.33 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: mmuCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.24 ","Gly  ","GGA  ","0.26 ","Gly  ","GGT  ","0.18 ","Gly  ","GGC  ","0.33 ", &
"Glu  ","GAG  ","0.60 ","Glu  ","GAA  ","0.40 ","Asp  ","GAT  ","0.44 ","Asp  ","GAC  ","0.56 ", &
"Val  ","GTG  ","0.46 ","Val  ","GTA  ","0.12 ","Val  ","GTT  ","0.17 ","Val  ","GTC  ","0.25 ", &
"Ala  ","GCG  ","0.10 ","Ala  ","GCA  ","0.23 ","Ala  ","GCT  ","0.29 ","Ala  ","GCC  ","0.38 ", &
"Arg  ","AGG  ","0.22 ","Arg  ","AGA  ","0.21 ","Ser  ","AGT  ","0.15 ","Ser  ","AGC  ","0.24 ", &
"Lys  ","AAG  ","0.61 ","Lys  ","AAA  ","0.39 ","Asn  ","AAT  ","0.43 ","Asn  ","AAC  ","0.57 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.16 ","Ile  ","ATT  ","0.34 ","Ile  ","ATC  ","0.50 ", &
"Thr  ","ACG  ","0.11 ","Thr  ","ACA  ","0.29 ","Thr  ","ACT  ","0.25 ","Thr  ","ACC  ","0.35 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.49 ","Cys  ","TGT  ","0.48 ","Cys  ","TGC  ","0.52 ", &
"End  ","TAG  ","0.23 ","End  ","TAA  ","0.28 ","Tyr  ","TAT  ","0.43 ","Tyr  ","TAC  ","0.57 ", &
"Leu  ","TTG  ","0.13 ","Leu  ","TTA  ","0.06 ","Phe  ","TTT  ","0.44 ","Phe  ","TTC  ","0.56 ", &
"Ser  ","TCG  ","0.05 ","Ser  ","TCA  ","0.14 ","Ser  ","TCT  ","0.20 ","Ser  ","TCC  ","0.22 ", &
"Arg  ","CGG  ","0.19 ","Arg  ","CGA  ","0.12 ","Arg  ","CGT  ","0.09 ","Arg  ","CGC  ","0.17 ", &
"Gln  ","CAG  ","0.75 ","Gln  ","CAA  ","0.25 ","His  ","CAT  ","0.40 ","His  ","CAC  ","0.60 ", &
"Leu  ","CTG  ","0.40 ","Leu  ","CTA  ","0.08 ","Leu  ","CTT  ","0.13 ","Leu  ","CTC  ","0.20 ", &
"Pro  ","CCG  ","0.10 ","Pro  ","CCA  ","0.28 ","Pro  ","CCT  ","0.31 ","Pro  ","CCC  ","0.31 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: rnoCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.24 ","Gly  ","GGA  ","0.25 ","Gly  ","GGT  ","0.17 ","Gly  ","GGC  ","0.34 ", &
"Glu  ","GAG  ","0.61 ","Glu  ","GAA  ","0.39 ","Asp  ","GAT  ","0.43 ","Asp  ","GAC  ","0.57 ", &
"Val  ","GTG  ","0.47 ","Val  ","GTA  ","0.11 ","Val  ","GTT  ","0.16 ","Val  ","GTC  ","0.25 ", &
"Ala  ","GCG  ","0.10 ","Ala  ","GCA  ","0.22 ","Ala  ","GCT  ","0.28 ","Ala  ","GCC  ","0.39 ", &
"Arg  ","AGG  ","0.21 ","Arg  ","AGA  ","0.20 ","Ser  ","AGT  ","0.15 ","Ser  ","AGC  ","0.24 ", &
"Lys  ","AAG  ","0.62 ","Lys  ","AAA  ","0.38 ","Asn  ","AAT  ","0.41 ","Asn  ","AAC  ","0.59 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.15 ","Ile  ","ATT  ","0.33 ","Ile  ","ATC  ","0.52 ", &
"Thr  ","ACG  ","0.11 ","Thr  ","ACA  ","0.28 ","Thr  ","ACT  ","0.24 ","Thr  ","ACC  ","0.37 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.50 ","Cys  ","TGT  ","0.45 ","Cys  ","TGC  ","0.55 ", &
"End  ","TAG  ","0.22 ","End  ","TAA  ","0.28 ","Tyr  ","TAT  ","0.40 ","Tyr  ","TAC  ","0.60 ", &
"Leu  ","TTG  ","0.13 ","Leu  ","TTA  ","0.06 ","Phe  ","TTT  ","0.42 ","Phe  ","TTC  ","0.58 ", &
"Ser  ","TCG  ","0.06 ","Ser  ","TCA  ","0.14 ","Ser  ","TCT  ","0.19 ","Ser  ","TCC  ","0.23 ", &
"Arg  ","CGG  ","0.20 ","Arg  ","CGA  ","0.12 ","Arg  ","CGT  ","0.09 ","Arg  ","CGC  ","0.18 ", &
"Gln  ","CAG  ","0.75 ","Gln  ","CAA  ","0.25 ","His  ","CAT  ","0.39 ","His  ","CAC  ","0.61 ", &
"Leu  ","CTG  ","0.41 ","Leu  ","CTA  ","0.08 ","Leu  ","CTT  ","0.12 ","Leu  ","CTC  ","0.20 ", &
"Pro  ","CCG  ","0.11 ","Pro  ","CCA  ","0.28 ","Pro  ","CCT  ","0.30 ","Pro  ","CCC  ","0.31 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: sceCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.12 ","Gly  ","GGA  ","0.22 ","Gly  ","GGT  ","0.47 ","Gly  ","GGC  ","0.19 ", &
"Glu  ","GAG  ","0.30 ","Glu  ","GAA  ","0.70 ","Asp  ","GAT  ","0.65 ","Asp  ","GAC  ","0.35 ", &
"Val  ","GTG  ","0.19 ","Val  ","GTA  ","0.21 ","Val  ","GTT  ","0.39 ","Val  ","GTC  ","0.21 ", &
"Ala  ","GCG  ","0.11 ","Ala  ","GCA  ","0.29 ","Ala  ","GCT  ","0.38 ","Ala  ","GCC  ","0.22 ", &
"Arg  ","AGG  ","0.21 ","Arg  ","AGA  ","0.48 ","Ser  ","AGT  ","0.16 ","Ser  ","AGC  ","0.11 ", &
"Lys  ","AAG  ","0.42 ","Lys  ","AAA  ","0.58 ","Asn  ","AAT  ","0.59 ","Asn  ","AAC  ","0.41 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.27 ","Ile  ","ATT  ","0.46 ","Ile  ","ATC  ","0.26 ", &
"Thr  ","ACG  ","0.14 ","Thr  ","ACA  ","0.30 ","Thr  ","ACT  ","0.35 ","Thr  ","ACC  ","0.22 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.30 ","Cys  ","TGT  ","0.63 ","Cys  ","TGC  ","0.37 ", &
"End  ","TAG  ","0.23 ","End  ","TAA  ","0.47 ","Tyr  ","TAT  ","0.56 ","Tyr  ","TAC  ","0.44 ", &
"Leu  ","TTG  ","0.29 ","Leu  ","TTA  ","0.28 ","Phe  ","TTT  ","0.59 ","Phe  ","TTC  ","0.41 ", &
"Ser  ","TCG  ","0.10 ","Ser  ","TCA  ","0.21 ","Ser  ","TCT  ","0.26 ","Ser  ","TCC  ","0.16 ", &
"Arg  ","CGG  ","0.04 ","Arg  ","CGA  ","0.07 ","Arg  ","CGT  ","0.15 ","Arg  ","CGC  ","0.06 ", &
"Gln  ","CAG  ","0.31 ","Gln  ","CAA  ","0.69 ","His  ","CAT  ","0.64 ","His  ","CAC  ","0.36 ", &
"Leu  ","CTG  ","0.11 ","Leu  ","CTA  ","0.14 ","Leu  ","CTT  ","0.13 ","Leu  ","CTC  ","0.06 ", &
"Pro  ","CCG  ","0.12 ","Pro  ","CCA  ","0.41 ","Pro  ","CCT  ","0.31 ","Pro  ","CCC  ","0.16 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: xlaCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.21 ","Gly  ","GGA  ","0.35 ","Gly  ","GGT  ","0.21 ","Gly  ","GGC  ","0.23 ", &
"Glu  ","GAG  ","0.48 ","Glu  ","GAA  ","0.52 ","Asp  ","GAT  ","0.57 ","Asp  ","GAC  ","0.43 ", &
"Val  ","GTG  ","0.36 ","Val  ","GTA  ","0.17 ","Val  ","GTT  ","0.27 ","Val  ","GTC  ","0.20 ", &
"Ala  ","GCG  ","0.07 ","Ala  ","GCA  ","0.32 ","Ala  ","GCT  ","0.33 ","Ala  ","GCC  ","0.27 ", &
"Arg  ","AGG  ","0.22 ","Arg  ","AGA  ","0.28 ","Ser  ","AGT  ","0.18 ","Ser  ","AGC  ","0.20 ", &
"Lys  ","AAG  ","0.49 ","Lys  ","AAA  ","0.51 ","Asn  ","AAT  ","0.52 ","Asn  ","AAC  ","0.48 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.23 ","Ile  ","ATT  ","0.42 ","Ile  ","ATC  ","0.35 ", &
"Thr  ","ACG  ","0.09 ","Thr  ","ACA  ","0.35 ","Thr  ","ACT  ","0.30 ","Thr  ","ACC  ","0.26 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.39 ","Cys  ","TGT  ","0.50 ","Cys  ","TGC  ","0.50 ", &
"End  ","TAG  ","0.18 ","End  ","TAA  ","0.43 ","Tyr  ","TAT  ","0.51 ","Tyr  ","TAC  ","0.49 ", &
"Leu  ","TTG  ","0.16 ","Leu  ","TTA  ","0.11 ","Phe  ","TTT  ","0.56 ","Phe  ","TTC  ","0.44 ", &
"Ser  ","TCG  ","0.05 ","Ser  ","TCA  ","0.16 ","Ser  ","TCT  ","0.23 ","Ser  ","TCC  ","0.19 ", &
"Arg  ","CGG  ","0.12 ","Arg  ","CGA  ","0.12 ","Arg  ","CGT  ","0.12 ","Arg  ","CGC  ","0.13 ", &
"Gln  ","CAG  ","0.64 ","Gln  ","CAA  ","0.36 ","His  ","CAT  ","0.50 ","His  ","CAC  ","0.50 ", &
"Leu  ","CTG  ","0.30 ","Leu  ","CTA  ","0.10 ","Leu  ","CTT  ","0.19 ","Leu  ","CTC  ","0.14 ", &
"Pro  ","CCG  ","0.09 ","Pro  ","CCA  ","0.37 ","Pro  ","CCT  ","0.32 ","Pro  ","CCC  ","0.22 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: ecoCFT = &
    RESHAPE( (/&
"Gly  ","GGG  ","0.16 ","Gly  ","GGA  ","0.15 ","Gly  ","GGT  ","0.34 ","Gly  ","GGC  ","0.35 ", &
"Glu  ","GAG  ","0.33 ","Glu  ","GAA  ","0.67 ","Asp  ","GAT  ","0.64 ","Asp  ","GAC  ","0.36 ", &
"Val  ","GTG  ","0.34 ","Val  ","GTA  ","0.17 ","Val  ","GTT  ","0.29 ","Val  ","GTC  ","0.20 ", &
"Ala  ","GCG  ","0.31 ","Ala  ","GCA  ","0.24 ","Ala  ","GCT  ","0.19 ","Ala  ","GCC  ","0.26 ", &
"Arg  ","AGG  ","0.05 ","Arg  ","AGA  ","0.08 ","Ser  ","AGT  ","0.17 ","Ser  ","AGC  ","0.23 ", &
"Lys  ","AAG  ","0.27 ","Lys  ","AAA  ","0.73 ","Asn  ","AAT  ","0.52 ","Asn  ","AAC  ","0.48 ", &
"Met  ","ATG  ","1.00 ","Ile  ","ATA  ","0.14 ","Ile  ","ATT  ","0.49 ","Ile  ","ATC  ","0.37 ", &
"Thr  ","ACG  ","0.24 ","Thr  ","ACA  ","0.19 ","Thr  ","ACT  ","0.19 ","Thr  ","ACC  ","0.38 ", &
"Trp  ","TGG  ","1.00 ","End  ","TGA  ","0.31 ","Cys  ","TGT  ","0.47 ","Cys  ","TGC  ","0.53 ", &
"End  ","TAG  ","0.09 ","End  ","TAA  ","0.60 ","Tyr  ","TAT  ","0.60 ","Tyr  ","TAC  ","0.40 ", &
"Leu  ","TTG  ","0.13 ","Leu  ","TTA  ","0.15 ","Phe  ","TTT  ","0.59 ","Phe  ","TTC  ","0.41 ", &
"Ser  ","TCG  ","0.13 ","Ser  ","TCA  ","0.15 ","Ser  ","TCT  ","0.17 ","Ser  ","TCC  ","0.14 ", &
"Arg  ","CGG  ","0.12 ","Arg  ","CGA  ","0.07 ","Arg  ","CGT  ","0.34 ","Arg  ","CGC  ","0.34 ", &
"Gln  ","CAG  ","0.66 ","Gln  ","CAA  ","0.34 ","His  ","CAT  ","0.59 ","His  ","CAC  ","0.41 ", &
"Leu  ","CTG  ","0.46 ","Leu  ","CTA  ","0.04 ","Leu  ","CTT  ","0.12 ","Leu  ","CTC  ","0.10 ", &
"Pro  ","CCG  ","0.47 ","Pro  ","CCA  ","0.21 ","Pro  ","CCT  ","0.19 ","Pro  ","CCC  ","0.14 "/), (/3,64/) )

  CHARACTER(LEN=5),DIMENSION(3,64) :: ppaCFT = &
    RESHAPE( (/&
"Phe  ","TTT  ","0.54 ","Phe  ","TTC  ","0.46 ","Ser  ","TCT  ","0.29 ","Ser  ","TCC  ","0.20 ", &
"Ser  ","TCA  ","0.19 ","Ser  ","TCG  ","0.09 ","Ser  ","AGT  ","0.15 ","Ser  ","AGC  ","0.09 ", &
"Tyr  ","TAT  ","0.46 ","Tyr  ","TAC  ","0.55 ","Cys  ","TGT  ","0.65 ","Cys  ","TGC  ","0.35 ", &
"Leu  ","TTA  ","0.16 ","Leu  ","TTG  ","0.33 ","Leu  ","CTT  ","0.16 ","Leu  ","CTC  ","0.08 ", &
"Leu  ","CTA  ","0.11 ","Leu  ","CTG  ","0.16 ","End  ","TAA  ","0.53 ","End  ","TGA  ","0.18 ", &
"End  ","TAG  ","0.29 ","Trp  ","TGG  ","1.00 ","Pro  ","CCT  ","0.35 ","Pro  ","CCC  ","0.15 ", &
"Pro  ","CCA  ","0.41 ","Pro  ","CCG  ","0.09 ","His  ","CAT  ","0.57 ","His  ","CAC  ","0.43 ", &
"Arg  ","CGT  ","0.16 ","Arg  ","CGC  ","0.05 ","Arg  ","CGA  ","0.10 ","Arg  ","CGG  ","0.05 ", &
"Arg  ","AGA  ","0.48 ","Arg  ","AGG  ","0.16 ","Gln  ","CAA  ","0.61 ","Gln  ","CAG  ","0.39 ", &
"Ile  ","ATT  ","0.50 ","Ile  ","ATC  ","0.30 ","Ile  ","ATA  ","0.19 ","Thr  ","ACT  ","0.40 ", &
"Thr  ","ACC  ","0.25 ","Thr  ","ACA  ","0.24 ","Thr  ","ACG  ","0.11 ","Asn  ","AAC  ","0.51 ", &
"Asn  ","AAT  ","0.49 ","Lys  ","AAA  ","0.47 ","Lys  ","AAG  ","0.53 ","Met  ","ATG  ","1.00 ", &
"Val  ","GTT  ","0.42 ","Val  ","GTC  ","0.23 ","Val  ","GTA  ","0.15 ","Val  ","GTG  ","0.19 ", &
"Ala  ","GCT  ","0.45 ","Ala  ","GCC  ","0.26 ","Ala  ","GCA  ","0.23 ","Ala  ","GCG  ","0.06 ", &
"Asp  ","GAT  ","0.58 ","Asp  ","GAC  ","0.42 ","Gly  ","GGT  ","0.44 ","Gly  ","GGC  ","0.14 ", &
"Gly  ","GGA  ","0.32 ","Gly  ","GGG  ","0.10 ","Glu  ","GAA  ","0.57 ","Glu  ","GAG  ","0.43 "/), (/3,64/) )

! SOLUTIONS

  TYPE DNA
    CHARACTER(LEN=9999) :: DNAseq=''  ! the actual DNA sequence,in ACGT nts
    INTEGER :: NumOlaps=0 
!    INTEGER :: NumOlaps=0              ! the total number of overlaps
    INTEGER :: OlapsPos(999,2)         ! the positions of the first and last
                                       !   nucleotides in the overlap
    INTEGER(KIND=1) :: NUMseq(9999)           ! the nt sequence as numbers (-1,-3,3,1)
    INTEGER(KIND=1) :: prot2cod(3333)          ! PROTpos to codon (1-64)
    INTEGER(KIND=1) :: nt2cod(9999)           ! DNApos to codon (1-64) or 0
    REAL :: MeltT(999)                 ! melting temps for the overlaps
    REAL :: TScore(999)                ! overlap-based score of MeltTm deviance
    REAL :: CScore(3333)               ! codon-based score of codon frequency
    REAL :: LScore(9999)              ! nt-based score of oligo length
    INTEGER :: RScore(9999)           ! nt-based score of repeats
    INTEGER :: PScore(9999)           ! nt-based score of pattern matching
    INTEGER :: MScore(9999)           ! nt-based score of potential mispriming
    INTEGER :: AScore(9999)           ! nt-based score of AT content
    INTEGER :: GScore(9999)           ! nt-based score of GC content
    INTEGER :: FScore(9999)           ! nt-based score of gap-fixed positions
    REAL :: TotalGScore=0              ! total score for GC content
    REAL :: TotalAScore=0              ! total score for AT content
    REAL :: TotalLScore=0              ! total score for oligo length
    REAL :: TotalCScore=0              ! total score for codons
    REAL :: TotalTScore=0              ! total score for temperature
    REAL :: TotalRScore=0              ! total score for repeats
    REAL :: TotalPScore=0              ! total score for patterns
    REAL :: TotalMScore=0              ! total score for mispriming
    REAL :: TotalFScore=0              ! total score for gap-fixed positions
    REAL :: OverallScore=0             ! Sum of all the total scores
    INTEGER :: RN=0             ! number of tandem repeats
    INTEGER :: RS1(9999)         ! starting position for primary seq
    INTEGER :: RS2(9999)         ! starting position for secondary seq
    INTEGER :: RLn(9999)          ! size of repeat (not oligo ends)
!    INTEGER(KIND=1) :: RX(9999)          ! direct=1,inverse=-1
    INTEGER :: RX(9999)          ! direct=1,inverse=-1
    INTEGER :: MN=0          ! number of potential misprimes
    INTEGER :: M1(9999)      ! starting position for potential misprime in prim
    INTEGER :: M2(9999)      ! starting position for potential misprime in seco
    INTEGER :: MX(9999)      ! Type of potential misprime (DS,IS,DA,IA)
    INTEGER :: MSN=0         ! number of actual misprimes
    INTEGER :: MS1(9999)     ! starting position for actual misprime in prim
    INTEGER :: MS2(9999)     ! starting position for actual misprime in seco
    INTEGER :: MSX(9999)     ! Type of actual misprime (DS,IS,DA,IA)
    INTEGER :: MOL(9999)     ! overlap the misprime is in
    INTEGER :: ntID_GC(9999)    ! window of GC content
    INTEGER :: ntID_AT(9999)    ! window of AT content
    INTEGER :: ntID_Tip(9999)   ! unique number for Tip matching
    INTEGER :: ntID_TipRC(9999) ! unique number for Tip (reverse complement)
    INTEGER :: ntID_Rep(9999)   ! unique number for Repeat matching
    INTEGER :: ntID_RepRC(9999) ! repeat matching (reverse complement)
    LOGICAL :: GapFixPos(9999)  ! should nt be fixed within a gap?
    LOGICAL :: Degen(9999)    ! true if the nt is degenerate
    INTEGER :: DegenNum(9999)   ! numerical index for degenerate sequence (1-11)
  END TYPE

  TYPE(DNA) :: CurrDNA
  TYPE(DNA) :: StoreDNA
  TYPE(DNA) :: BestDNA
  TYPE(DNA) :: BestOverlapDNA

END MODULE dnaworks_data
