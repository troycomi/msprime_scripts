#!/bin/bash

mod=$1
pop=$2

[[ -d version_1 ]] && rm -rf version_1

python ../Admixture_Simulation.py -p $pop \
    -m $mod --out-dir version_1 \
    -l $3 -s 2 -n 0.1 -d 0.1 \
    -g AF_EU_1.5e-5 -G AF_B_14e-5
