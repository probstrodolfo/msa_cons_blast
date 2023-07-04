# Rodolfo Probst v.1.1 June 19, 2023
# Licensed to PSF under a Contributor Agreement.

"""This Python script compiles consensus from original MSA. But it does more!
The structure is a bit messy, but idea is to transform FASTQ into FASTA;
Then filter data based on quality control and length (check this section!);
Then align MSA (using MAFFT). Then generating consensus FASTA;
Then BLASTing consensus and parsing xml into txt with the top 10 sequences output of BLAST.

Make sure to check information below!

Script runs with the subprocess module for the MAFFT command with the --auto option. 
Change stringency as you pleased (gap penalty, etc.).
Make sure you have the MAFFT software installed on your system
and accessible in the command-line environment for this code to work correctly.

IMPORTANT!: Make sure to replace input file as your FASTQ (line  43)
"""

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from collections import Counter
import os
import subprocess


# Configure Entrez with your email address - required step
Entrez.email = "probstrodolfo@gmail.com"

# Select FASTQ file and filter based on length (blocked for now)
#input_fastq = "Formica_subelongata_EXHS0005.fastq"
#output_filtered_fasta = "filtered_sequences.fasta"
#min_length = 450
#max_length = 700

# Select FASTQ file, change Phred if needed (start with 10), make notes about # used
# Filter shorter sequences and change input_fastq to represent your filename

input_fastq = "Formica_subnitens_EXHS0027.fastq"
output_filtered_fasta = "filtered_sequences.fasta"
quality_threshold = 6
min_length = 300

# Load the FASTQ file and filter sequences based on quality score
filtered_records = []
for record in SeqIO.parse(input_fastq, "fastq"):
    if min(record.letter_annotations["phred_quality"]) >= quality_threshold and (len(record.seq) >= min_length):
        filtered_records.append(record)

# Save the filtered sequences as a new FASTA file
SeqIO.write(filtered_records, output_filtered_fasta, "fasta")

# Perform MAFFT alignment using the input FASTA MSA
input_fasta = output_filtered_fasta
mafft_output = os.path.splitext(input_fasta)[0] + "_alignment.fasta"

# Run MAFFT alignment with custom gap penalties
mafft_cmd = f"mafft --auto --op 5 --ep 1 {input_fasta} > {mafft_output}"
subprocess.call(mafft_cmd, shell=True)

# Load the multiple alignment
alignment = AlignIO.read(mafft_output, "fasta")

# Create an empty consensus sequence
consensus_seq = ""

# Iterate over each position in the alignment
for i in range(alignment.get_alignment_length()):
    # Count the occurrences of each nucleotide or amino acid
    column = alignment[:, i]
    counts = Counter(column)

    # Find the most frequent element(s) at this position
    max_count = max(counts.values())
    consensus_elements = [element for element, count in counts.items() if count == max_count]

    # Choose one consensus element (randomly or using additional criteria if needed)
    consensus_element = consensus_elements[0]

    # Append the consensus element to the consensus sequence
    consensus_seq += consensus_element

# Format the consensus sequence as a plain string
consensus_seq = str(consensus_seq)

# Create the consensus alignment file name
consensus_alignment_file = os.path.splitext(input_fasta)[0] + "_consensus.fasta"

# Save the consensus sequence as a new FASTA file
with open(consensus_alignment_file, "w") as file:
    file.write(f">{os.path.splitext(input_fastq)[0]}_consensus\n{consensus_seq}\n")

print(f"Consensus sequence saved as {consensus_alignment_file}")

# Perform BLAST search
result_handle = NCBIWWW.qblast("blastn", "nt", consensus_seq)

# Save BLAST results to an XML file
blast_xml_output_file = "blast_results.xml"
with open(blast_xml_output_file, "w") as blast_xml_file:
    blast_xml_file.write(result_handle.read())
    
print(f"BLAST results saved as {blast_xml_output_file}")
    
    
# Close the result handle to release the resources
result_handle.close()

# Parse the BLAST results from the XML file
blast_records = NCBIXML.parse(open(blast_xml_output_file))

# Save top 10 BLAST results to a text file
blast_txt_output_file = "blast_results.txt"
with open(blast_txt_output_file, "w") as blast_txt_file:
    for record in blast_records:
        blast_txt_file.write(f">Sequence: {record.query}\n")
        count = 0
        for alignment in record.alignments:
            if count >= 10:
                break
            blast_txt_file.write(f"Hit: {alignment.title}\n")
            for hsp in alignment.hsps:
                blast_txt_file.write(f"Score: {hsp.score}\n")
                blast_txt_file.write(f"Alignment length: {hsp.align_length}\n")
                blast_txt_file.write(f"Query: {hsp.query}\n")
                blast_txt_file.write(f"Subject: {hsp.sbjct}\n")
            blast_txt_file.write("\n")
            count += 1

print(f"Top 10 BLAST results saved as {blast_txt_output_file}. You are done! Go enjoy your barcodes!")
