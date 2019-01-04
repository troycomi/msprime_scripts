#!/bin/bash

module load anaconda3
conda activate msprime_scripts

srun --ntasks=1 --cpus-per-task=10 --time=01:10:00 --mem=1G snakemake --use-conda -j 10

exit
snakemake --cluster-config 'della_cluster.yaml' \
    --cluster "sbatch --cpus-per-task={cluster.n} \
                --mem={cluster.memory} --time={cluster.time} \
                --output=slurm_out/%x_%A --job-name={cluster.jobname}" \
    --use-conda -w 60 -rp -j 250
