from Bio import SeqIO  # Import SeqIO module from the Bio package
import sys  # Import sys module for command-line arguments

# Get the input FASTA file name and output BED file name from command-line arguments
fasta_input_name = sys.argv[1]
bed_output_name = sys.argv[2]

# Open the input FASTA file in "read" mode and output BED file in "write" mode
with open(fasta_input_name, "r") as fasta_file:
    with open(bed_output_name, "w") as bed_file:
        current_chromosome = None  # Variable to keep track of the current chromosome
        for record in SeqIO.parse(fasta_file, "fasta"):
            chromosome_name = record.id  # Retrieve the ID (name) of the chromosome from the FASTA record
            if chromosome_name != current_chromosome:
                current_chromosome = chromosome_name
                start_position = 0  # Start position of the current chromosome
            for position in range(len(record.seq)):
                end_position = position + 1  # Calculate the end position based on the current position
                nucleotide = record.seq[position]  # Retrieve the nucleotide at the current position
                # Write the chromosome name, start position, end position, and nucleotide to the BED file
                bed_file.write("{}\t{}\t{}\t{}\n".format(chromosome_name, start_position, end_position, nucleotide))
                start_position = end_position  # Update the start position for the next iteration
