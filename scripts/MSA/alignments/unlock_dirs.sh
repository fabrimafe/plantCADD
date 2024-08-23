name_file=$1

module load snakemake

while IFS= read -r line; do
    cd ${line}_pipeline
    snakemake --unlock
    cd ..
done < ${name_file}