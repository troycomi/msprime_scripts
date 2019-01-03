#!/bin/bash

module load anaconda3
conda activate msprime_scripts

snakemake --cluster-config 'della_cluster.yaml' \
    --cluster "sbatch --cpus-per-task={cluster.n} \
                --mem={cluster.memory} --time={cluster.time} \
                --output=slurm_out/%x_%A --job-name={cluster.jobname}" \
    --use-conda -w 60 -rp -j 250
