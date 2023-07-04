# msa_cons_blast

-----------------------------------------------------------------------------
msa_cons_blast 1.0 Readme - Jun/19/2023

Rodolfo Probst (probstrodolfo@gmail.com)

(a) Science Research Initiative (SRI) - College of Science, University of Utah
(b) School of Biological Sciences, University of Utah

-----------------------------------------------------------------------------
This Python program compiles consensus from original MSA. But it does more!
The structure is a bit messy, but idea is to transform FASTQ file into FASTA;
Then filter data based on quality control and length (check this section in code);
Then align MSA (using MAFFT by penalizing gappy regions). Then generating consensus FASTA;
Then BLASTing consensus and parsing xml into txt with the top 10 sequences output of BLAST.

Make sure to check information below!

Script runs with the subprocess module for the MAFFT command with the --auto option. 
Change stringency as you please (gap penalty, etc.).
Make sure you have the MAFFT software installed on your system
and accessible in the command-line environment for this code to work correctly (python version >3.0)
