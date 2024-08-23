#!/bin/bash
module load tabix
module load miniconda/4.10.3_environmentally
conda activate NN_new

# parsing of raw files
bash 1_parse_repeats.sh

# Go over all chromosomes
for chr in {0..12}; do

    # Expand to bed
    bash 2_expand_to_bed.sh ${chr}

    # Dummy up
    python dummy_beds.py ${chr}

    # Add header hashtag
    sed -i '1s/^/#/' repeats_chr_${chr}_dummied.bed

    # Gzip
    mv repeats_chr_${chr}_dummied.bed repeats_feature_chr_${chr}.bed
    bgzip repeats_feature_chr_${chr}.bed
   
    # Remove tmps
    rm parsed_repeats_chr_${chr}.mod repeats_chr_${chr}.bed
done