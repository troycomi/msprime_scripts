#!/bin/bash

module load anaconda3
conda activate msprime_scripts

[[ -d version_1 ]] && rm -rf version_1

python ../Admixture_Simulation.py -p nonAfr \
    -m Tenn --out-dir version_1 \
    -l 1e4 -s 2 -n 0.1 -d 0.1 > teststd.out

exit
./comp_samples.sh
