#!/bin/bash
# A script to generate codon positions
#module load miniconda/22.11.1_environmentally
conda activate NN_new

gff_file=$1
chromosomes=$2

# Make the all codon position file to concatenate into
rm -f codon_pos/codon_pos_all_chrs.bed
touch codon_pos/codon_pos_all_chrs.bed
for i in {1..3}
do
    # Generate the specific codon position
    #python scripts/fab_codon_pos.py ${gff_file} "codon_pos/tmp_${gff_file}_codon_${i}.bed" --pos ${i}
    # Add the codon value
    awk -v codon_val="$i" '{print $0 "\t" codon_val}' codon_pos/tmp_${gff_file}_codon_${i}.bed > codon_pos/${gff_file}_codon_${i}.bed
    # Concatenate them all in order to sepearte in chromosomes
    cat "codon_pos/${gff_file}_codon_${i}.bed" >> codon_pos/codon_pos_all_chrs.bed
    rm "codon_pos/${gff_file}_codon_${i}.bed"
done

cd codon_pos
# Split it into chromosomes
awk '{print > $1}' codon_pos_all_chrs.bed
# Change their names
for i in ${chromosomes}
do
    sort -k1,1n -k2,2n ${i} > cds_codon_pos_chr_${i}.bed
done
cd ..
rm codon_pos/codon_pos_all_chrs.bed