#!/bin/bash

module load hdf5/gcc/1.8.16
module load anaconda
conda activate msprime_scripts2

PYFI=msprime_scripts/msprime.Admixture_Simulate.py

[[ -s test.out ]] && rm test.out 
python $PYFI -p nonAfr -o Tenn \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c haplo &>>test.out
python $PYFI -p nonAfr -o Tenn  \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c F4Dstat &>>test.out
python $PYFI -p nonAfr -o Tenn  \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c vcf &>>test.out
python $PYFI -p nonAfr -o Sriram \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c haplo &>>test.out
python $PYFI -p nonAfr -o SplitPop \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c haplo &>>test.out
python $PYFI -p nonAfr -o test \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c haplo &>>test.out
python $PYFI -p nonAfr -o Tenn_nomod \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c haplo &>>test.out
python $PYFI -p nonAfr -o Tenn_pulsed \
    -l 1e4 -s 2 -i 2 -n 0.1 -d 0.1 -c haplo &>>test.out

echo DONE
