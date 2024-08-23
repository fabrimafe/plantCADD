from Bio import SeqIO
import sys
import os

# Get args
fasta_input_name = sys.argv[1]
output_directory = sys.argv[2]

# Open the fasta file to read from
with open(fasta_input_name, "r") as fasta_file:
    # Go over fasta chromosomes
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Get chromosome name
        chromosome_name = record.id
        # Create output file path
        output_filename = os.path.join(output_directory, f"fasta_nuc_CAP_chr_{chromosome_name}.bed")

        # Open the bed file to write to
        with open(output_filename, "w") as bed_file:
            # Define current chromosome
            current_chromosome = None
            # Initialize start position
            start_position = 0

            # Write out the chromosome into the bed file
            for position in range(len(record.seq)):
                # Offset for bed format
                end_position = position + 1
                # Get the nucleotide for that position
                nucleotide = record.seq[position].upper()  # Capitalize the nucleotide
                # Write into the output file
                bed_file.write("{}\t{}\t{}\t{}\n".format(chromosome_name, start_position, end_position, nucleotide))
                # Advance start position
                start_position = end_position
