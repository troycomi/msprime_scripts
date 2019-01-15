#!/bin/bash

module load anaconda3
conda activate msprime_scripts

THREADS=8
SNAKE_DIR=/tigress/tcomi/abwolf_abc/.snakemake

#to view dependency graph (use only a few seeds!)
#snakemake --dag | dot -Tsvg > dag.svg

#download all environments needed. 
#Note: shadow prefix determines location of snakemake
#data, including environments, must match srun below
snakemake --use-conda --create-envs-only \
    --shadow-prefix $SNAKE_DIR 

#unlock if job was interrupted
#snakemake --unlock

#perform workflow.  
srun --ntasks=1 --cpus-per-task=$THREADS --time=01:00:00 \
    --mem=4G snakemake --use-conda -j $THREADS \
    --shadow-prefix  $SNAKE_DIR
    #--rerun-incomplete
