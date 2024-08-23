# This script is able to handle gz bed files as well.

import sys
import math
import gzip
import subprocess

def open_bed_file(bed_file):
    # Open the BED file
    if bed_file.endswith('.gz'):
        # If the BED file is gzipped, use gzip module to open it
        return gzip.open(bed_file, 'rt')
    else:
        # Otherwise, open it as a regular text file
        return open(bed_file, 'r')

def parse_gff(gff_file, chromosome):
    # List to store gene positions for the chromosome
    gene_positions = []
    
    # Read the GFF file to retrieve gene positions for the chromosome
    with open(gff_file, 'r') as gff:
        for line in gff:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                # Check if the feature is a CDS and matches the chromosome
                if len(fields) > 2 and fields[2] == 'gene' and fields[0] == chromosome:
                    gene_start = int(fields[3])
                    gene_end = int(fields[4])
                    strand = fields[6]
                    gene_positions.append((gene_start, gene_end, strand))
    return gene_positions

def calculate_distances(gff_file, bed_file, chromosome):
    # Parse the GFF file to retrieve gene positions for the chromosome
    gene_positions = parse_gff(gff_file, chromosome)
    
    # Print the header
    print(f"#chr\tstart\tend\tupstream\tdownstream")
    
    # Open the BED file for reading
    with open_bed_file(bed_file) as bed:
        # Iterate over each line in the BED file
        for line in bed:
            fields = line.strip().split('\t')
            if fields[0] == chromosome:
                position = int(fields[1])
                closest_upstream = math.inf
                closest_downstream = math.inf
                for gene_start, gene_end, strand in gene_positions:
                    # Calculate distances based on strand
                    if strand == '+':
                        if gene_end < position and (position - gene_end) < closest_downstream:
                            closest_downstream = position - gene_end
                        elif gene_start > position and (gene_start - position) < closest_upstream:
                            closest_upstream = gene_start - position
                        elif gene_start <= position <= gene_end:
                            closest_upstream = closest_downstream = 0
                            break
                    elif strand == '-':
                        if gene_start > position and (gene_start - position) < closest_upstream:
                            closest_upstream = gene_start - position
                        elif gene_end < position and (position - gene_end) < closest_downstream:
                            closest_downstream = position - gene_end
                        elif gene_start <= position <= gene_end:
                            closest_upstream = closest_downstream = 0
                            break
                    
                # Print the results
                print(f"{chromosome}\t{position}\t{position+1}\t{closest_upstream}\t{closest_downstream}")

# Check if the number of command-line arguments is correct
if len(sys.argv) != 4:
    print("Usage: python script.py input.gff input.bed chromosome")
    sys.exit(1)

gff_file = sys.argv[1]
bed_file = sys.argv[2]
chromosome = sys.argv[3]

# Calculate distances between positions and the closest genes for the specified chromosome
calculate_distances(gff_file, bed_file, chromosome)
