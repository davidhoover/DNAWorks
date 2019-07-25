DNAWorks
========

Automatic oligonucleotide design for PCR-based gene synthesis 

 DNAWorks v3.2.4
 David Hoover
 May 04, 2017
 
 DNAWorks takes as input nucleotide and/or protein sequences, codon
 information, and other variables, and attempts to optimize a synthetic
 gene.  It then outputs the gene with a variety of histograms and metrics
 for judging the probability of success for generating the gene by PCR.  It
 also outputs the oligonucleotide sequences required for PCR synthesis of
 the synthetic gene.
 
 This program is based on this publication:
 
   Hoover DM, Lubkowski J. DNAWorks: an automated method for designing
   oligonucleotides for PCR-based gene synthesis. Nucleic Acids Res. 2002 May
   15;30(10):e43. PubMed PMID: 12000848; PubMed Central PMCID: PMC115297.
 
 Kindly reference this publication if you use this for your work.
 

Installation
============

Currently, DNAWorks is written in Fortran.  It will require a Fortran compiler on a UNIX system.

If you do not have gfortran, make, or git, then on a Linux machine, install these packages (Ubuntu):

```
apt-get install gfortran make git_hub
```

or on Centos

```
yum install gfortran make git
```

Then download DNAWorks and compile with make:

```
git clone https://github.com/davidhoover/DNAWorks.git
cd DNAWorks
make
```
 
and the dnaworks executable should compile.

Run
===

Instructions (can be displayed by typing ```./dnaworks -help```):

```
 COMMAND-LINE OPTIONS
 ==============================================================================
 
 The command line is as follows:
 
   % dnaworks [ inputfile ] [ -t0 | -t1 | -t2 | -t3 ]
 
 The default inputfile is 'DNAWORKS.inp'.  All options, except for those
 on the command line, are read from the inputfile.  See below for a complete
 description of the options.
 
 The flags -t0, -t1, -t2, and -t3 are for testing purposes.  They report
 the internal actions within the program based on the level input.
 
   -t0   Relatively simple output, only subroutine names
   -t1   Most subroutine names reported
   -t2   Heavy output, all subroutines, some functions
   -t3   Way too much output, all subroutines and functions reported
 
 INPUTFILE OPTIONS
 ==============================================================================
 
 The input is case insensitive, except for quoted strings.  Any string
 can be quoted, but it's not necessary unless the case must be preserved or
 if there are spaces or special characters (#,!).  The quotes can
 be single or double, but must begin and end around the intended
 string.  
 
 Any text that follows a '#' or '!' is considered comments, and will
 be ignored.
 
 Options in the inputfile are of the following types:
 
   [ S ]   string
 
 Strings are converted to uppercase, unless quoted (either " or '')
 
   [ #I ]  integer number
   [ #R ]  real number
 
 Integers are, well, integers.  Real numbers can be floating point numbers
 (e.g., 12.345) or scientific notation (e.g., -12.36E+4).
 
   [ name ]  directive
 
 Directives are special strings the enable or disable particular functions.
 In general, only the first 4 or 5 characters are actually read, so they
 can be abbreviated.
 
 Directives must be placed flat against the left margin of the input file,
 otherwise they will be ignored.
 
 ------------------------------------------------------------------------------
 
 INPUT DIRECTIVES:
 
   [ tbio ]
 
 The method of gene synthesis employed by DNAWorks is termed 
 'thermodynamically balanced', in that all the oligonucleotides should 
 assemble and anneal at the same temperature.  The amplification occurs 
 everywhere at once, and ideally can generate the gene with just one round 
 of PCR.  However, there are sticky cases where the gene does not amplify, 
 and constructing the gene in pieces is not successful.
 
 A more controlled method of gene synthesis, termed 'thermodynamically 
 balanced inside-out', was developed for cases where problems occurred 
 during PCR synthesis (Gao, et al., 2003). In an assembly set of 
 oligonucleotides, the first half of the oligos are all synthesized in the 
 sense orientation, and the other half are synthesized as reverse complements 
 in the anti-sense orientation of the gene. The gene assembly and amplification 
 is thus done in steps of 0.4-0.6 kb from the center pair of 
 oligonucleotides outward.
 
 Enabling tbio will enable thermodynamically balanced inside-out output.
 
 
   [ nogaps ]
 
 
 By default, DNAWorks will try to keep all oligos the same size as the chosen
 length.  If the size is beyond the sizes required for the chosen Tm, gaps
 are introduced between overlap regions.  The directive nogaps will keep oligos
 as short as possible, with no gaps between the overlap regions.
 
 Restricting oligos to no gaps may slow down the optimization somewhat, and 
 may result in higher scores due to a higher probability of misprimes.
 
 ------------------------------------------------------------------------------
 
 INPUT OPTIONS:
 
 
   logfile [ S ]
 
 
 The default output file is 'LOGFILE.txt'.  Entering a string after the
 logfile option will change the name of the logfile.
  
  
   title [ S ]
 
 
 It's always good to give the output a title to keep it unique and to give
 you an easy way to keep track of what the output is.
 
 
   timelimit [ #I ]
 
 
 Set a time limit for the run, in seconds.  This keeps the program from 
 running forever.  A value of 0 (the default) means no limit.
 
 
   solutions [ #I ]
 
 
 Normally DNAWorks only generates a single solution for a set of parameters.
 Since the optimization involves a lot of random number calls, and that it is
 impossible to get to the 'true minimum' by Monte Carlo methods, sometimes
 generating more than one solutions is a good thing.  Look for the best
 solution in the end.  The range is 1-99.
 
 
   melting [ #I ] [ low #I high #I ] [ tolerance #I ]
 
 
 This governs the chosen melting or annealing temperature for the oligos.
 Giving a single integer (between 55 and 75) will generate a single solution.
 A range of melting temperatures can be given with the low and high options,
 and a solution for each temperature will be generated.  The tolerance value
 is by default +/- 1 degree, but it can be modified.  Don't set it too high
 or the point of the program can be lost!
 
 
   length [ #I ] [ low #I high #I ] [ random ]
 
 
 This sets the ideal length of the oligo.  Because the oligos can have gaps,
 they can be as long as you wish, but remember that errors accumulate in
 synthetic DNA oligos very quickly beyond around 50 nts!
 
 By default, an attempt is made to force all oligos to be the same size as the
 chosen length.  On occasion this can lead to a higher probability of 
 misprimes.  Also, this can limit successful optimization when sequences 
 are gapfixed (see below), since gap position and size will be limited.  In
 this case, enabling the length directive random causes oligos to be 
 designed with random length (between 20 nt and the length chosen).
 
 
   frequency [ threshold #I ] [ random ] [ strict ] [ score ]
 
 
 The frequency threshold is the cutoff for which codons are used for
 reverse translation of protein sequences into DNA.  For example, a value of
 20 will allow only those codons whose frequencies equal or exceed 20%.
 
 By default, DNAWorks uses the highest frequency codons for the initial
 reverse translation of the protein sequences.  Having the random option
 present causes the program to choose the initial codons at random.
 
 By default, DNAWorks always uses the two highest frequency codons for 
 optimization.  To override this default, enabling strict will 
 force the program to strictly use only those codons that are within the 
 chosen codon frequency threshold.  Be careful, because setting a high 
 codon frequency threshold (>20%) and strict will result in many protein 
 residues with a single codon available, and thus very little room for 
 optimization.
 
 To accelerate convergence, DNAWorks does not continuously score codon 
 frequency. This is allowed because only the highest frequency codons are 
 usually used.  However, for the particularly picky user, enabling scored will 
 force the program to continuously evaluate the codon frequency score. This 
 will have the effect of increasing the overall frequency of codons (at 
 the cost of other scores...). 
 
 
   concentration [ oligo #R ] [ sodium #R ] [ magnesium #R ]
 
 
 The concentration of oligonucleotides, monovalent cations (Na+, K+), and 
 magnesium in the PCR reaction can have profound effects on the annealing 
 temperatures of the oligonucleotides.  The user can enter the desired 
 concentrations for the PCR reaction.
 
 The effects of these components on the annealing temperature is based on 
 the program HyTher (Nicolas Peyret, Pirro Saro and John SantaLucia, Jr.).
 
 Values are in moles per liter, and can be entered in scientific notation 
 for simplicity.
 
 Oligonucleotides must be between 100 um (1E-4 M) and 1 nm (1E-9 M), 
 monovalent cations must be between 10 and 1000 mM, and magnesium must be 
 between 0 and 200 mM.
 
 
   repeat [ #I ]
 
 
 DNAWorks continuously monitors the synthetic gene for any repeats that 
 occur within the gene.  A repeat can be a direct repeat, an inverted
 repeat (which can result in a hairpin), or a palindromic repeat.  If a 
 repeat occurs that is above a certain length, it can lead to stable 
 annealing of oligos to unexpected positions and mispriming.  Such mispriming
 can result in either no PCR product, or a long smear on a gel.
 
 The value for repeat governs the minimum length of nucleotides considered
 a repeat.  The default value is 8.  Increasing this number will
 decrease the number of repeats found, while decreasing it will do the
 opposite.
 
 
   misprime [ #I ] [ tip #I ] [ max #I ]
 
 
 The major flaw to PCR-based gene synthesis is mispriming.  This occurs when
 an oligo anneals to an unexpected position on the PCR template.  To prevent
 this from happening, DNAWorks compares the ends of each oligo with the
 current synthetic sequence and analyzes its potential to anneal to that
 site.  
 
 A misprime is a special variant of a repeat, in that it only occurs at the
 business end (3') of an oligo.  
 
 The first number for misprime is the length of the sequence to compare.  The 
 default value is 18.
 
 The tip number is number of nucleotides that must be exactly identical at
 the tip of the oligo.  The default is 6.  This value is based on little more
 than guessing, but increasing it will cause very few misprimes to be
 identified, and decreasing will cause too many to be identified.
 
 The max number is the maximum number of non-identical nucleotides in the
 misprime sequence.  The default is 8.  This number is again a guess.  It
 is generally not understood why non-identical sequences anneal to each other,
 but it is based on structural and electrostatic principles that are way too
 difficult to incorporate into this program.  Again, increasing the number
 results in too many misprimes to be identified, decreasing it causes too few.
 
 Needless to say, the misprime value is just plain prudence, but not 
 necessarily fact.
 
 
   weight [ twt #R ] [ cwt #R ] [ rwt #R ] [ mwt #R ] [ gwt #R ] [ awt #R ]
     [ lwt #R ] [ pwt #R ] [ fwt #R ]
 
 
 DNAWorks optimizes a synthetic gene by evaluating the scores of a set of 
 features: annealing temperature (T), codon frequency (C), repeat (R), 
 misprime potential (M), GC- (G) and AT- (A) content, length (L), gapfix (F)
 and pattern constraining (P).  The default weights of each individual feature 
 score are set to 1.  By increasing the weight of an individual feature, the 
 final output can be nudged to favoring one feature over the others.  For 
 example, in the case where the potential synthetic genes for a set of 
 sequences chronically suffers from high number of repeats, increasing the 
 weight of the repeat score (RWT) might decrease the final repeat score at 
 the expense of the other feature scores.
 
 Beware, as modulating the weights is not fully tested.  Remember that this 
 merely skews the results toward one feature or another, and may do more 
 harm than good.  In most cases keeping the weights balanced is the best 
 approach.
 
 
   previous [ #I ] [ S ]
 
 
 DNAWorks allows old sets of oligonucleotides to be read back with a new,
 mutant gene.  It then calculates scores for the mutant gene with overlap
 positions and parameters identical to the original solution.  It then
 outputs only those oligonucleotides that need to be changed.  This is very
 useful for generating mutants, since in general only one or two new oligos
 need to be synthesized.
 
 The integer refers to the previous solution number, and the string is the
 name of the previous logfile.
 
 ------------------------------------------------------------------------------
   
 INPUT SECTIONS:
 
   nucleotide [ reverse | gapfix ]
     ...
   //
 
 Nucleotide sequences can only include A,C,G, or T in the nucleotide
 section.  They can also include degenerate sequences:
 
     B = C or G or T         rev. compl. = V
     D = A or G or T         rev. compl. = H
     H = A or C or T         rev. compl. = D
     K = G or T              rev. compl. = M
     M = A or C              rev. compl. = K
     N = A or C or G or T    rev. compl. = N
     R = A or G              rev. compl. = Y
     S = C or G              rev. compl. = S
     V = A or C or G         rev. compl. = B
     W = A or T              rev. compl. = W
     Y = C or T              rev. compl. = R
 
   protein [ reverse | gapfix ]
     ...
   //
 
 Protein sequences can be input through the protein section, but can only 
 include the single-letter abbreviations of the 20 standard amino acids 
 (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y).  Stop codons are designated by X.
 
 The reverse directive causes the nucleotide sequence (either original or
 translated from the protein sequence) to be reversed on incorporation in
 the synthetic gene.
 
 The gapfix directive is used when the sequence should not fall within 
 overlap regions, but rather only in the gaps or overhangs that are single
 stranded in the annealed assembly prior to PCR.  This is advantageous for
 subsequent mutations by oligonucleotide replacement. For example, if a 
 synthetic gene will be exhaustively mutated at a single codon, having the
 codon entirely within a gap region will allow its mutation by replacing a
 single oligonucleotide, rather than two or three.
 
 The gapfix directive will enable Fixed Gap Scoring.  Any nt that are
 designated as gapfixed but fall within overlap regions will increase the
 global score.  DNAWorks will then try to minimize the score by moving the
 gap regions toward the gapfixed nucleotides.  Because gap regions are 
 generally short (less than 10 nt), the sequence should be very short.  
 Otherwise the global score will remain quite high, and other features (Tm,
 repeats, misprimes) will not receive as much attention.
 
 Gapfixing is much more effective when oligo lengths are allowed be 
 randomized, rather than fixed to the length chosen by default.  See
 length option, above, for more details.
 
 
   codon [ ecoli2 | E. coli | C. elegans | D. melanogaster | H. sapiens | 
     M. musculus | R. novegicus | S. cerevesiae | X. laevis | P. pastoris ]
     [ ...
   // ]
 
 Codon frequencies can be entered manually in the codon section using 
 GCG-format codon frequencies.  If a directive corresponding to a given 
 organism is present, the codon frequency for that organism will be used.
 
   pattern
     ...
   //
 
 Nucleotide patterns can be screened if entered in the pattern section.
 Pattern sequences can be normal or degenerate nucleotide sequences.
```


Good luck!
