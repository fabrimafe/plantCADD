#!/bin/bash
module load bedtools
gff_file="/home/labs/alevy/omerbar/backups/S_lycopersicum_data/SL5/SL5.gff3"
scaler="/home/labs/alevy/omerbar/NN_general/OLD_GENE_DIST_EFFED_dataset_M44_standardized_full_dataset_old_filt_13/4.17M_standard_scaler.joblib"
num_of_muts=$1
# Concatenate all mini prediction files
if [ -e "ALL_PREDICTIONS.bed" ]; then\
    echo "Continue..."
else
    echo "Concatenate all predictions..."
    cat PREDICTIONS_chr_1.bed >> ALL_PREDICTIONS.bed
    for i in {2..12}
    do
        tail -n +2 PREDICTIONS_chr_${i}.bed >> ALL_PREDICTIONS.bed
    done
fi

echo "Get positions of top 1%"
# Get the top 1% of the predictions
python top_delet_get_pos.py ALL_PREDICTIONS.bed top_1_percent_deleterious.bed ${num_of_muts}

# Break into chromosomes
awk '{print > "top_1_percent_deleterious_chr_"$1".bed"}' top_1_percent_deleterious.bed

echo "Get features for positions"
for i in {1..12}
do
    echo "Processing chr ${i}"
    # Intersect with features (randomly chosen feature file - prediction file from this set)
    bedtools intersect -a 28_test_model_loss_0.5657972918848411.pth_predictions_chr_${i}.txt -b top_1_percent_deleterious_chr_${i}.bed -wa -sorted> top_1_percent_deleterious_chr_${i}_features.bed

    # Retain wanted features
    python extract_features.py top_1_percent_deleterious_chr_${i}_features.bed ${scaler} top_1_percent_deleterious_chr_${i}_features_extracted.bed 
done

# Put all together
echo "Concatenate all features"
cat top_1_percent_deleterious_chr_1_features_extracted.bed.csv >> ALL_DELET_POS_FEATURES.bed
for i in {2..12}
do
    tail -n +2 top_1_percent_deleterious_chr_${i}_features_extracted.bed.csv >> ALL_DELET_POS_FEATURES.bed
done

echo "Get transcript ids for positions"
# Get 0based positions for mrna from gff
cat  ${gff_file}  | \
awk '{if($3=="mRNA"){print $1"\t"$4-1"\t"$5"\t"$9"\t"$5-$4-1}}' | grep -v "^#" | \
awk '{if($1!=0)print}' | \
awk -F'\t' '{gsub("ID=", "", $4);split($4,a,";");print $1"\t"$2"\t"$3"\t"a[1]"\t"$5}' | \
sort -k1,1n -k2,2n > mrna_positions_gff.bed

# Get number of cols for feature file
num_cols=$(head -n 1 ALL_DELET_POS_FEATURES.bed | awk '{print NF}')

# Add header to the file and add transcript_id
printf "%s\ttranscript_id\n" "$(head -1 ALL_DELET_POS_FEATURES.bed)" > ALL_DELET_POS_FEATURES_MRNA.bed
printf "%s\ttranscript_id\n" "$(head -1 ALL_DELET_POS_FEATURES.bed)" > ALL_DELET_POS_FEATURES_REST_GENOME.bed

echo "Get gene ids for positions"
# intersect to get the gene ids feature into the main table
bedtools intersect -a ALL_DELET_POS_FEATURES.bed -b mrna_positions_gff.bed -wa -wb -sorted | \
awk -v cols=${num_cols} '{for(i=1;i<=cols;i++) printf "%s\t", $i; printf "%s\n", $(NF-1)}' >> ALL_DELET_POS_FEATURES_MRNA.bed
# For each position without a gene, use the min dist from up/down and get a position that corresponds to a gene, add that

# Subtract to get the rest of the deleterious variants

bedtools subtract -a ALL_DELET_POS_FEATURES.bed -b mrna_positions_gff.bed| \
awk -v cols=${num_cols} '{for(i=1;i<=cols;i++) printf "%s\t", $i; printf ".\n"}' >> ALL_DELET_POS_FEATURES_REST_GENOME.bed

# Get closest gene to the rest of the positions
printf "%s\tgff_chr\tgff_start\tgff_end\tgff_transcript_id\ttranscript_length\n" "$(head -1 ALL_DELET_POS_FEATURES_REST_GENOME.bed)" > TMP_ALL_DELET_POS_FEATURES_REST_GENOME_T_ID.bed
bedtools closest -a ALL_DELET_POS_FEATURES_REST_GENOME.bed -b mrna_positions_gff.bed >> TMP_ALL_DELET_POS_FEATURES_REST_GENOME_T_ID.bed

# Remove transcript duplications
python remove_dups.py TMP_ALL_DELET_POS_FEATURES_REST_GENOME_T_ID.bed ALL_DELET_POS_FEATURES_REST_GENOME_T_ID.bed

# Combine all
cat ALL_DELET_POS_FEATURES_MRNA.bed  <(tail -n +2 ALL_DELET_POS_FEATURES_REST_GENOME_T_ID.bed) | sort -k1,1n -k2,2n > FINAL_DELET_POS.bed

rm top_1_percent_deleterious.bed ALL_DELET_POS_FEATURES.bed ALL_DELET_POS_FEATURES_MRNA.bed ALL_PREDICTIONS.bed
rm ALL_DELET_POS_FEATURES_REST_GENOME.bed mrna_positions_gff.bed top_1_percent_deleterious_chr_#chr.bed  top_1_percent_deleterious_chr_1_features_extracted.bed.csv
rm TMP_ALL_DELET_POS_FEATURES_REST_GENOME_T_ID.bed ALL_DELET_POS_FEATURES_REST_GENOME_T_ID.bed
for i in {1..12};do rm top_1_percent_deleterious_chr_${i}_features_extracted.bed.csv top_1_percent_deleterious_chr_${i}.bed top_1_percent_deleterious_chr_${i}_features.bed;done