#!/bin/bash

#module load anaconda3
conda activate msprime_scripts

set -euo pipefail

##Modify to change number of jobs to be run at once
THREADS=8

##to view dependency graph (use only a few seeds!)
#snakemake --dag | dot -Tsvg > dag.svg

##download all environments needed. 
snakemake --use-conda --create-envs-only

##unlock if job was interrupted
#snakemake --unlock

##perform workflow.  
srun --ntasks=1 --cpus-per-task=$THREADS --time=02:00:00 \
    --mem=4G snakemake --use-conda -j $THREADS \
    -pr \
    #-R $(snakemake --list-input-changes) \
    #--rerun-incomplete
