#!/bin/bash
# A script to generate codon position in a chromosomal bed format
module load tabix
module load miniconda/22.11.1_environmentally
conda activate NN_new

gff_file=$1
org_name=${gff_file%.*}
# Per codon position possible
for i in {1..3}; do 
    # Get all nucleotides for this codon position
    python fab_codon_pos.py  ${gff_file} ${org_name}_nucs_codon_${i}.bed_TMP --pos ${i} 
    # Append the codon position value per file
    awk -v pos_val="$i" '{print $0"\t"pos_val}' ${org_name}_nucs_codon_${i}.bed_TMP > ${org_name}_nucs_codon_${i}.bed
    # Remove tmp files
    rm ${org_name}_nucs_codon_${i}.bed_TMP
done

echo "Combining all position files..."
# Combine all positions 
rm -f ${org_name}_all_codon_positions.bed_TMP; touch ${org_name}_all_codon_positions.bed_TMP
for i in {1..3}; do
    cat ${org_name}_nucs_codon_${i}.bed >> ${org_name}_all_codon_positions.bed_TMP
    # Remove tmp files
    rm ${org_name}_nucs_codon_${i}.bed
done

echo "Sorting..."
# Sort by chromosome/start
sort -n -k1,1n -k2,2n ${org_name}_all_codon_positions.bed_TMP > ${org_name}_all_codon_positions.bed

python take_min_from_dups.py ${org_name}_all_codon_positions.bed CLEAN_${org_name}_all_codon_positions.bed

# The point of the command below is to split the file based on the first column's values, the chromosome name, and also bgzip it, all in one go.
awk -v name=${org_name} '{filename = name "_codons_chr_" $1 ".bed"; print  | "bgzip > " filename ".gz"}' CLEAN_${org_name}_all_codon_positions.bed

# Remove tmp files
rm ${org_name}_all_codon_positions.bed_TMP