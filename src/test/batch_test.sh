#!/bin/bash

module load anaconda3
conda activate msprime_scripts

set -e

for len in 1e5 1e6
do
    for mod in Tenn Sriram SplitPop Tenn_nomod Tenn_pulsed
    do
        for pop in EAS EUR nonAfr AFR modHum
        do
            if [[ ! -d version_2/${len}_${pop}_${mod} ]]; then
                echo using model $mod with population $pop
                ./get_samples.sh $mod $pop $len

                ./comp_samples.sh $mod $pop $len
                echo no diffs...
            fi
        done
    done
done
echo done

