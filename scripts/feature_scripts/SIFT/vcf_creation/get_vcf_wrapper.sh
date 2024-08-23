#!/bin/bash
# Create the chromosomal vcfs

# Args
gff=$1
fasta=$2
output_name=$3 # .vcf will be added
mem=$4

mkdir -p bsub_logs
# Get chromosome list
chromosomes=$( cat ${fasta} | grep -oE '^>[^[:space:]]+' | sed 's/^>//' )
for i in ${chromosomes}
do
    bsub -n 1 -m public_hosts -q new-long -J sl5_vcf_sift_generation_${i} -e sl5_vcf_sift_generation_${i}_%J.error -o sl5_vcf_sift_generation_${i}_%J.output -R rusage[mem=${mem}] "python create_vcf.py ${gff} ${i} ${fasta} ${output_name}_chr_${i}.vcf"
done
#---------------------------------------------------------------
# PT 2 to be done after PT1 is done
# Concatenate them all#
#rm -f ${output_name}_all.vcf
#touch ${output_name}_all.vcf
#for i in ${chromosomes}
#do

#    cat ${output_name}_chr_${i}.vcf | sort -k1,1n -k2,2n >> ${output_name}_all.vcf
#    rm ${output_name}_chr_${i}.vcf
#done

#mv *.output *.error bsub_logs/