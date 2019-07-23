#!/bin/bash

module load anaconda3
conda activate msprime_scripts

set -euo pipefail

# Sample script to run multiple iterations of snakemake
for N1 in 0.3 0.4
do
    # make config file to overwrite default config
    ### NOTE ###
    # due to limits in --config flag, can't change nested values
    # due to limits in --configfile flag, can't use multiple files
    echo '---
paths:
    admixed_dir: "__BASE_OUTPUT__/n1_'${N1}'_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0"

msprime:
    n1: '$N1'' > tempconfig.yaml

    ## perform workflow with sbatch
    snakemake --cluster-config della_cluster.yaml \
        --cluster "sbatch --cpus-per-task={cluster.n} \
            --mem={cluster.memory} --time={cluster.time} \
            --output=slurm_out/%x_%A --job-name={cluster.jobname} \
            --parsable -A eeb" \
        --use-conda \
        -pr \
        -w 60 -j 50 \
        --configfile tempconfig.yaml
done
