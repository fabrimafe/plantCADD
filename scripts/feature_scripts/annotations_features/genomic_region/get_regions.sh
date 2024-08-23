#!/bin/bash
# Generate genomic regions as dummies, not regarding silly mRNA overlaps.
# This way we will parse through hirarchy of genomic regions later on, when we collapse this feature for plotting.
module load bedtools
module load miniconda
conda activate NN_new
#This has to be run before the script, because it is divided by chromosomes################
chromosome_sizes="/home/labs/alevy/omerbar/features/expanded_SL5/feature_file_generation/chromosome.sizes"
gff_file="/home/labs/alevy/omerbar/features/expanded_SL5/gff_features/SL5.gff3"
awk '{print > "SL5_gff_chr_"$1".gff3"}' ${gff_file}
awk '{print > "chr_"$1".size"}' ${chromosome_sizes}

chr=$1
# Expands the gff to bed, makes it 0 based, and adds the genomic regions as dummy features
python parse_gff_to_bed.py SL5_gff_chr_${chr}.gff3 SL5_regions_chr_${chr}.bed

# Get intergenic regions
bedtools complement -i SL5_gff_chr_${chr}.bed -g SL5_gff_chr_${chr}.bed > SL5_gff_chr_${chr}_intergenic.bed

# Expand the regions
fasta_file="/home/labs/alevy/omerbar/features/expanded_SL5/feature_file_generation/SL5_chr_${chr}.bed.gz"
bedtools intersect -a <(zcat ${fasta_file}) -b  SL5_gff_chr_${chr}_intergenic.bed -wa -wb -sorted | awk '{print $1"\t"$2"\t"$3"\t"1}' > intergenic_regions_chr_${chr}.bed
bgzip intergenic_regions_chr_${chr}.bed
rm SL5_gff_chr_${chr}_intergenic.bed chr_${chr}.size SL5_gff_chr_${chr}.gff3

