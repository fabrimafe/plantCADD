#!/bin/bash
# Intersecting feature file with novel and neutral mutations, picking the correct lines from CDS regions (to match the mutation)
module load bedtools
module load tabix

# Flags

novel=false
neutral=false

while getopts ":ab" opt; do
  case ${opt} in
    a)
      novel=true
      ;;
    b)
      neutral=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))


chr=$1
novel_dir=$2
neut_dir=$3
feature_file=/home/labs/alevy/omerbar/features/arabidopsis/main_table/AT_features_SIFT_chr_${chr}.bed.gz
novel_file=./${novel_dir}novel_mut_chr_${chr}.bed
neut_file=./${neut_dir}neut_mut_chr_${chr}.bed
num_cols_feat=$(awk 'NR==2 {print NF; exit}' <(zcat ${feature_file}))
num_cols_nov_mut=$(awk 'NR==2 {print NF; exit}' ${novel_file})
num_cols_neut_mut=$(awk 'NR==2 {print NF; exit}' ${neut_file})

echo $neutral 
echo $novel
if [ "$neutral" = true ]; then
    echo "collect neutral"
    # grep -v '-' to remove deletions 
    # awk if for the cutoff of mappability quality

    bedtools intersect -a $feature_file -b ${neut_file} -wa -wb | \

    # For the neutral mutation dataset, we pick oppositely. We match based on the ancestral allele, and not the derived one (lycopersicum allele for neutral) 
    # $15 - genomic region is CDS
    # $5 - SIFT ALT allele (the mutation it predicts for)
    # $61 - anc_allele (always last col), the mutation identity for neutral mutations (intersect adds chr/start/stop of other file)
    # Zero here "0\n" is neutral label 0
    # We need to get the mutation allele because if there's no SIFT info, we won't have it. So we add the last col of the mutation file as well
    # Second awk is to remove duplicate lines. It goes by column values beside the SIFT's (last 6), so change accordingly
    # Last print in here is the alt allele, useful when it is not CDS and we do not have that info from SIFT
    awk -v feat_cols=$num_cols_feat -v neut_cols=$num_cols_neut_mut '{if (($15 == 1 && $5 == $61) || ($15 == 0)) {for (i = 1; i <= feat_cols; i++) printf "%s\t", $i; printf "%s\t0\n", $(feat_cols + neut_cols);}}' | \
    #awk '{line = $0; $NF=$(NF-1)=$(NF-2)=$(NF-3)=$(NF-4)=$(NF-5)=""; if (!seen[$0]++) {print line}}' | \
    sed 's/\r//g' | bgzip > ./${neut_dir}neut_features_chr_${chr}.bed.gz
fi

if [ "$novel" = true ]; then
    echo "collect novel"
    # First awk - filter out coding regions that do not match the novel mutation (each cds position has 4 rows, collapsed into one)
    # Second awk -  remove duplicate lines. It goes by column values beside the SIFT's (last 6), so change accordingly
    bedtools intersect -a <(zcat $feature_file) -b ${novel_file} -wa -wb -sorted | \
    # Only get rows where the mutation nucleotide is the same as the ALT allele (features mutation_allele and ALT_ALLELE)
    # Last print in here is the alt allele, useful when it is not CDS and we do not have that info from SIFT
    awk -v feat_cols=$num_cols_feat -v nov_cols=$num_cols_nov_mut '{if (($15 == 1 && $5 == $61) || ($15 == 0)) {for (i = 1; i <= feat_cols; i++) printf "%s\t", $i; printf "%s\t1\n", $(feat_cols + nov_cols);}}' | \
    #awk '{line = $0; $NF=$(NF-1)=$(NF-2)=$(NF-3)=$(NF-4)=$(NF-5)=""; if (!seen[$0]++) {print line}}' | \
    # Removing the fucking bullshit ^M dos line end signatures
    sed 's/\r//g' | bgzip > ./${novel_dir}novel_features_chr_${chr}.bed.gz
fi
