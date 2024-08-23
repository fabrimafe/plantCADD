run_all_chrs_gene_dist.sh is a wrapper to get all the chromosome's gene distance.
It generates upstream and downstream gene distance based on the closest coordinate of a gene considering the strands(whether the start/end for whatever the strand is)
The second argument for the wrapper is the skeleton feature file (chr/start/end/ref_nucleotide)

ex.
bash generate_gene_dist.sh /home/labs/alevy/omerbar/backups/S_lycopersicum_data/SL5/SL5.gff3 /home/labs/alevy/omerbar/features/expanded_SL5/feature_file_generation/SL5_chr_11.bed.gz 11