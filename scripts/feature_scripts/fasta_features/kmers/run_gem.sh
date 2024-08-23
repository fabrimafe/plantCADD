#!/bin/bash
# The point is to get a position wise precentage of unique kmers per all the kmers available in which the position is in the window of. 

# Modules
module load bedtools # this is for the expansion
module load bedops # this is for the wig2bed
module load tabix

# Args
myfasta=$1 #data/Arabidopsis_thaliana.fa
myoutput=$2 #mappability/map
kmer=$3 #35

# Getting chromosome names
chromosomes=$( cat ${myfasta} | grep -oE '^>[^[:space:]]+' | sed 's/^>//' )
# Getting fasta file organism name
fasta_base_name=$(basename "$myfasta")
fasta_base_name="${fasta_base_name%%.*}"

# Index
echo "Indexing..."
gem-indexer -i ${myfasta} -o ${myoutput} # Output gets a prefix .gem
# Mappability
echo "Generating mappability..."
gem-mappability -I ${myoutput}.gem -l ${kmer} -o ${myoutput}${kmer} # Output gets a prefix .mappability
# Gem to wig
echo "Converting mappability to wig..."
gem-2-wig -I ${myoutput}.gem -i ${myoutput}${kmer}.mappability -o ${myoutput}
# Wig to bed
echo "Converting wig to bed..."
wig2bed < "${myoutput}.wig" > ${myoutput}.bed
# Expand the bed file after converting it from wig
echo "Parsing the bed file..."
awk '{for(i=$2;i<$3;i++) print $1"\t"i"\t"i+1"\t"$5}' ${myoutput}.bed > ${myoutput}_all_chr_${kmer}kmers.bed
# Expand the bed file to chromosomes
awk '{print > $1".bed"}' ${myoutput}_all_chr_${kmer}kmers.bed

mkdir -p kmer_results
# rename the files
for i in ${chromosomes}
do
	mv ${i}.bed ${fasta_base_name}_chr_${i}_${kmer}kmers.bed
	bgzip ${fasta_base_name}_chr_${i}_${kmer}kmers.bed
	mv ${fasta_base_name}_chr_${i}_${kmer}kmers.bed kmer_results/
done

echo "Done!"
# Remove useless stuff?
rm "${myoutput}.log" "${myoutput}.gem" "${myoutput}.wig" "${myoutput}${kmer}.mappability" "${myoutput}.bed" "${myoutput}.sizes"