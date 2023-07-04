# msa_cons_blast

-----------------------------------------------------------------------------
msa_cons_blast 1.0 Readme - Jun/19/2023

Rodolfo Probst (probstrodolfo@gmail.com)

(a) Science Research Initiative (SRI) - College of Science, University of Utah
(b) School of Biological Sciences, University of Utah

-----------------------------------------------------------------------------
This Python program compiles consensus from original MSA. But it does more!
The structure is a bit messy, but idea is to transform FASTQ file into FASTA;
THEN filter data based on quality control and length (check this section in code);
THEN align MSA (using MAFFT by penalizing gappy regions). Finally, after generating consensus FASTA,
consensus sequence will be BLASTed and xml saved and subsequently parsed into txt with the top 10 sequences output of BLAST query (based on % similarity).

Make sure to check information below!

. Script runs with the subprocess module for the MAFFT command with the --auto option. 

. Change stringency as you please (gap penalty, etc.)

. Make sure you have the MAFFT software installed on your system

. and accessible in the command-line environment for this code to work correctly (python version >3.0)

. You need to have Anaconda or Miniconda installed to use Bio

. For current version, you need to replace input fastq every time you run program (I am still working on parallelizing, but issue with unique quality scores for different MSAs impedes that for now).
