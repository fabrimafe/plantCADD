#!/bin/bash
# A simple script to run all jobs easily
for i in 15 20 25 30 35
do
	cd /home/labs/alevy/omerbar/features/SL5/fasta_features/kmers/kmers${i}
	bsub -n 1 -m public_hosts -q new-long -J SL5_kmers${i} -e SL5_kmers${i}_%J.error -o SL5_kmers$i_%J.output -R rusage[mem=20000] "bash run_gem.sh /home/labs/alevy/omerbar/features/SL5/feature_file_generation/SL5.fasta SL5 ${i}"
done
