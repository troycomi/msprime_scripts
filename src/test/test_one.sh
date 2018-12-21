#!/bin/bash

module load anaconda3
conda activate msprime_scripts

set -e
len=1e3

mod=$1
pop=$2
if [[ $# -eq 2 ]] ; then
    echo using model $mod with population $pop
    ./get_samples.sh $mod $pop $len &
    cd ../../../old_v/msprime_scripts
    ./get_samples.sh $mod $pop $len &
    cd - > /dev/null
    wait
    ./comp_samples.sh $mod $pop
    echo no diffs...
else
    echo getting debug
    python ../Admixture_Simulation.py -p $pop \
        -m $mod --out-dir version_1 \
        -l $len -s 2 -n 0.1 -d 0.1 \
        -g AF_EU_1.5e-5 -G AF_B_14e-5 \
        --debug debug.txt --options options.txt
    cd ../../../old_v/msprime_scripts/output
    python ../msprime.Admixture_Simulate.py -p $pop \
        --migration_AF_EU 1.5e-5 --migration_AF_B 14e-5 \
        -o $mod -l $len -s 2 -n 0.1 -d 0.1 -c debug > debug.out 2> options.out
    cd - > /dev/null
fi

echo done

