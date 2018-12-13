#!/bin/bash

module load anaconda3
conda activate msprime_scripts

[[ -s test.out ]] && rm test.out Tenn*
python ../Admixture_Simulation.py -p nonAfr -o Tenn \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c haplo &>>test.out
python ../Admixture_Simulation.py -p nonAfr -o Tenn  \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c F4Dstat &>>test.out

./comp_samples.sh
