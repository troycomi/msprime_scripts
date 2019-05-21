#!/bin/bash

module load anaconda3
conda activate msprime_scripts
PATH=$PATH:/usr/licensed/anaconda3/2019.3/bin

set -euo pipefail

# make slurm_out if not exists
[[ -d slurm_out ]] || mkdir slurm_out

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
    --jobscript jobtemplate.sh \
    --configfile della_config.yaml
    #-R $(snakemake --list-input-changes) \
    #--rerun-incomplete
