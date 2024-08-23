#!/bin/bash
# The main bash script to run all the other subscripts
# THE FIRST SCRIPT ASSUMES THE FASTA AND ANC (BGZIPPED) FILES ARE IN THE MAIN FOLDER
# THE SECOND SCRIPT ASSUMES THE ANC STATE HAS 4 COLS AND SO DOES THE FASTA NUC
# THE THIRD SCRIPT ASSUMES THE VCF FILES ARE IN THE ./vcfs FOLDER. THIER NAMES SHOULD BE 'CHR_NAME.vcf.gz'
# THE FOURTH SCRIPT ASSUMES A GFF FILE WITH CDS POSITIONS IN THE MAIN FOLDER

# For quiet bsub commands for this script
export BSUB_QUIET=Y

# Modules
module load bedtools
conda activate NN_new

# Parameters
ref_genome=$1
anc_state=$2
gff_file=$3
home_dir=$(pwd)

# Making sure subdirectories exist
mkdir -p fasta_nucleotides
mkdir -p anc_states
mkdir -p anc_fasta_intersect
mkdir -p completely_fixed
mkdir -p nearly_fixed
mkdir -p rare_mutations
mkdir -p codon_pos
mkdir -p enrichment

# ------------------------------ANALYSIS TO DIFFERENT GROUPS---------------------
# These are done here because there's no point in scripting them apart, the analysis is different
chromosomes=$( cat ${ref_genome} | grep -oE '^>[^[:space:]]+' | sed 's/^>//' )

#                            --------------------Completely fixed--------------------
echo "Parsing completely fixed mutations..."
rm -f completely_fixed/completely_fixed_all_chrs_CODONS.bed
touch completely_fixed/completely_fixed_all_chrs_CODONS.bed
rm -f completely_fixed/ALL_completely_fixed_mutations.bed
touch completely_fixed/ALL_completely_fixed_mutations.bed
for i in ${chromosomes}
do  
    # This is the actual mutation bulk
    cat anc_fasta_intersect/anc_not_ref_chr_${i}.bed >> completely_fixed/ALL_completely_fixed_mutations.bed
    # This is for the enrichment check
    # Intersect anc!=ref files with codon position files and concatenate to master file
    # This group is not intersecting with vcfs as we need variants outside the vcfs
    bedtools intersect -a anc_fasta_intersect/anc_not_ref_chr_${i}.bed -b codon_pos/cds_codon_pos_chr_${i}.bed -sorted -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9}' >> completely_fixed/completely_fixed_all_chrs_CODONS.bed
done

#                            -----------------------------------------------------------
#                            -----------------------Nearly Fixed------------------------
echo "Parsing nearly fixed mutations..."
cd nearly_fixed
# Creating masters for both groups
rm -f nearly_fixed_anc_not_ref_5_all_chrs_CODONS.bed
touch nearly_fixed_anc_not_ref_5_all_chrs_CODONS.bed
rm -f nearly_fixed_anc_is_ref_95_all_chrs_CODONS.bed
touch nearly_fixed_anc_is_ref_95_all_chrs_CODONS.bed
rm -f ALL_nearly_fixed_mutations.bed
touch ALL_nearly_fixed_mutations.bed

for i in ${chromosomes}
do
    # First is intersecting with VCF to get pool of positions (bulk mutations)
    # Second is to intersect with codon position (enrichment check)
    # A. Anc==Ref, AF>95% in vcf
    bedtools intersect -a ../anc_fasta_intersect/anc_is_ref_chr_${i}.bed -b ../vcfs/chr_${i}.vcf.bed -sorted -wa -wb | awk '{if($11>0.95) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}'  > anc_is_ref_vcf_95_chr_${i}.bed
    cat anc_is_ref_vcf_95_chr_${i}.bed >> ALL_nearly_fixed_mutations.bed
    bedtools intersect -a anc_is_ref_vcf_95_chr_${i}.bed -b ../codon_pos/cds_codon_pos_chr_${i}.bed -sorted -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10}' >> nearly_fixed_anc_is_ref_95_all_chrs.bed

    # B. Anc!=Ref, AF<5% in vcf
    bedtools intersect -a ../anc_fasta_intersect/anc_not_ref_chr_${i}.bed -b ../vcfs/chr_${i}.vcf.bed -sorted -wa -wb | awk '{if($11<0.05) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}' > anc_not_ref_vcf_5_chr_${i}.bed
    cat anc_not_ref_vcf_5_chr_${i}.bed >> ALL_nearly_fixed_mutations.bed
    bedtools intersect -a anc_not_ref_vcf_5_chr_${i}.bed -b ../codon_pos/cds_codon_pos_chr_${i}.bed -sorted -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10}' >> nearly_fixed_anc_not_ref_5_all_chrs.bed
done
# Sorting the bulk mutation file
sort -k1,1n -k2,2n ALL_nearly_fixed_mutations.bed > ALL_nearly_fixed_mutations_sorted.bed
rm ALL_nearly_fixed_mutations.bed; mv ALL_nearly_fixed_mutations_sorted.bed ALL_nearly_fixed_mutations.bed
cd ..
#                            -----------------------------------------------------------
#                            ---------------------------Rare----------------------------
echo "Parsing rare mutations..."
cd rare_mutations
# Creating masters for both groups
rm -f rare_mutations_anc_not_ref_99_all_chrs.bed
touch rare_mutations_anc_not_ref_99_all_chrs.bed
rm -f rare_mutations_anc_is_ref_1_all_chrs.bed
touch rare_mutations_anc_is_ref_1_all_chrs.bed
rm -f ALL_rare_mutations.bed
touch ALL_rare_mutations.bed

for i in ${chromosomes}
do
    # First is intersecting with VCF to get pool of positions (bulk mutations)
    # Second is to intersect with codon position  (enrichment check)

    # A. Anc==Ref, AF<1% in vcf
    bedtools intersect -a ../anc_fasta_intersect/anc_is_ref_chr_${i}.bed -b ../vcfs/chr_${i}.vcf.bed -sorted -wa -wb | awk '{if($11<0.01) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}'  > anc_is_ref_vcf_1_chr_${i}.bed
    cat anc_is_ref_vcf_1_chr_${i}.bed >> ALL_rare_mutations.bed
    bedtools intersect -a anc_is_ref_vcf_1_chr_${i}.bed -b ../codon_pos/cds_codon_pos_chr_${i}.bed -sorted -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10}' >> rare_mutations_anc_is_ref_1_all_chrs.bed


    # B. Anc!=Ref, AF>99% in vcf
    bedtools intersect -a ../anc_fasta_intersect/anc_not_ref_chr_${i}.bed -b ../vcfs/chr_${i}.vcf.bed -sorted -wa -wb | awk '{if($11>0.99) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}' > anc_not_ref_vcf_99_chr_${i}.bed
    cat anc_not_ref_vcf_99_chr_${i}.bed >> ALL_rare_mutations.bed
    bedtools intersect -a anc_not_ref_vcf_99_chr_${i}.bed -b ../codon_pos/cds_codon_pos_chr_${i}.bed -sorted -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10}' >> rare_mutations_anc_not_ref_99_all_chrs.bed
done

# Sorting bulk mutation file
sort -k1,1n -k2,2n ALL_rare_mutations.bed > ALL_rare_mutations_sorted.bed
rm ALL_rare_mutations.bed; mv ALL_rare_mutations_sorted.bed ALL_rare_mutations.bed
cd ..
#                            -----------------------------------------------------------
# ------------------------------------------PRINTING RESULTS-------------------------------------------

# Printing out the different groups and saving them as well
rm -f enrichment/enrichment.txt
touch enrichment/enrichment.txt
echo "Completely fixed mutations:" >> enrichment/enrichment.txt
python scripts/print_enrichment.py completely_fixed/completely_fixed_all_chrs_CODONS.bed | tee -a enrichment/enrichment.txt
echo "Nearly fixed mutations:" >> enrichment/enrichment.txt
python scripts/print_enrichment.py <(cat nearly_fixed/nearly_fixed_anc_not_ref_5_all_chrs.bed nearly_fixed/nearly_fixed_anc_not_ref_95_all_chrs.bed) | tee -a enrichment/enrichment.txt
echo "Rare mutations:" >> enrichment/enrichment.txt
python scripts/print_enrichment.py <(cat rare_mutations/rare_mutations_anc_is_ref_1_all_chrs.bed rare_mutations/rare_mutations_anc_is_ref_99_all_chrs.bed) | tee -a enrichment/enrichment.txt

# Resetting quiet bsub commands
export BSUB_QUIET=N