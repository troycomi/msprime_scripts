#!/bin/bash

module load anaconda3
conda activate msprime_scripts

set -euo pipefail

## to view dependency graph (use only a few seeds!)
#snakemake --dag | dot -Tsvg > dag.svg

## unlock if job was interrupted
#snakemake --unlock

## perform workflow with sbatch
snakemake --cluster-config della_cluster.yaml \
    --cluster "sbatch --cpus-per-task={cluster.n} \
        --mem={cluster.memory} --time={cluster.time} \
        --output=slurm_out/%x_%A --job-name={cluster.jobname} \
        --parsable -A eeb" \
    --use-conda \
    -pr \
    -w 60 -j 50 \
    --configfile della_config.yaml
    #-R $(snakemake --list-input-changes) \
    #--rerun-incomplete
