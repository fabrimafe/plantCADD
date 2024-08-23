# To do this I need to load some modules:
module load miniconda/22.11.1_environmentally
module load snakemake

conda activate msa
job_name=$1
mem=$2

# And run the first pipeline command:
snakemake --use-conda -j 999 --cluster "bsub -n 1 -m public_hosts -q new-all.q -J ${job_name}_mini -e ${job_name}_mini_%J.err -o ${job_name}_mini_%J.out -R rusage[mem=${mem}]" --configfile config/config.yaml -R call_conservation
