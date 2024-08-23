from Bio import SeqIO  # Import SeqIO module from the Bio package
import sys  # Import sys module for command-line arguments

# Get the input FASTA file name and output sizes file name from command-line arguments
fasta_input_name = sys.argv[1]

# Open the input FASTA file in "read" mode and output sizes file in "write" mode
with open(fasta_input_name, "r") as fasta_file:
    with open("chromosomes.sizes", "w") as sizes_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            chromosome_name = record.id  # Retrieve the ID (name) of the chromosome from the FASTA record
            chromosome_size = len(record.seq)  # Calculate the size of the chromosome sequence
            # Write the chromosome name and size to the sizes file
            sizes_file.write("{}\t{}\n".format(chromosome_name, chromosome_size))
