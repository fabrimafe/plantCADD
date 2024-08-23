This feature is a dummy one
It will make as many dummies as categories for positions in the gff annotations file (gene,mRNA,exon,CDS, etc..). It will also add a gene_count for that feature in the file.
Two files per chromosome, SL5_regions_chr_X.bed.gz for what is described above, and intergenic_regions_chr_X.bed.gz noting the intergenic regions
When you append this feature, to do it for both.

example
bash get_regions.sh 1