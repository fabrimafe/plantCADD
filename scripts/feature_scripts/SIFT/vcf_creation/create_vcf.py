from Bio import SeqIO
import sys

gff_file = sys.argv[1]
chromosome_name = sys.argv[2]
fasta_file = sys.argv[3]
output_file = sys.argv[4]

# Load the reference genome from the FASTA file
genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Open the GFF file and output file
with open(gff_file) as gff, open(output_file, "w") as vcf:
    # Write the VCF header
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for line in gff:
        # Skip comments and empty lines
        if line.startswith("#") or not line.strip():
            continue

        # Split the line into columns
        columns = line.strip().split("\t")

        # Extract the feature type and chromosome
        feature_type = columns[2]
        chromosome = columns[0]

        # Process only CDS features on the specified chromosome
        if feature_type == "CDS" and chromosome == chromosome_name:
            start = int(columns[3]) + 1  # Add 1 to convert to 1-based position
            end = int(columns[4])

            # Retrieve the reference sequence for the chromosome
            reference_sequence = str(genome[chromosome_name].seq)

            for position in range(start, end + 1):
                ref_allele = reference_sequence[position - 1]

                # Duplicate each row and add ALT alleles A, C, T, G
                for alt_allele in ["A", "C", "T", "G"]:
                    vcf.write("{chromosome}\t{position}\t.\t{ref}\t{alt}\t.\t.\t.\n".format(
                        chromosome=chromosome, position=position, ref=ref_allele, alt=alt_allele))
