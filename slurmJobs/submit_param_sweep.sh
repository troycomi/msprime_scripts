#!/bin/bash

set -euo pipefail

module load anaconda3

OUTPUT_DIR=~/projects/abwolf_abc/msprime_scripts/slurmJobs/slurm_out/output/
[[ -d $OUTPUT_DIR ]] || mkdir $OUTPUT_DIR

export BASE_OPTIONS="-p nonAfr -m Tenn -e 1006 -a 2050 -l 1e7 --out-dir $OUTPUT_DIR"
export SWEEP_OPTIONS=$(python ../src/Parameter_Sweeper.py \
        -p "n;0.2;:.2f" \
        -p "d;0.1;:.2f")
        #-p "n;0:0.01:0.05,0.07,0.1:0.05:0.2;:.2f" \
        #-p "d;0:0.01:0.05,0.07,0.1;:.2f")

export MSPRIME_LOC=~/projects/abwolf_abc/msprime_scripts/src/Admixture_Simulation.py
#Optional, unset to use random seed
export SEED=2
export BATCH_SIZE=1

readarray -t SWEEP_ARR <<<"$SWEEP_OPTIONS"
numjobs=$((${#SWEEP_ARR[@]} / $BATCH_SIZE))
#check if remainder is needed
if (( $numjobs * $BATCH_SIZE == ${#SWEEP_ARR[@]} )); then
    numjobs=$(($numjobs - 1))
fi

JOB_ID_SIM=$(sbatch -D slurm_out \
    --job-name=${USER}-msprime \
    --parsable \
    --array=0-$numjobs \
    msprime.slurm)

echo Submitted job $JOB_ID_SIM with ${#SWEEP_ARR[@]} simulations in $(($numjobs + 1)) batches
