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
 
Good luck!
