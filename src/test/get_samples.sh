#!/bin/bash

mod=$1
pop=$2

python ../Admixture_Simulation.py -p $pop \
    -m $mod --out-dir version_2/${3}_${pop}_${mod} \
    -l $3 -s 2 -n 0.1 -d 0.1 \
    -g AF_EU_1.5e-5 -G AF_B_14e-5
