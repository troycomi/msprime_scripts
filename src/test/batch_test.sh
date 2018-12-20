#!/bin/bash

module load anaconda3
conda activate msprime_scripts

set -e
len=1e5

for mod in Tenn Sriram SplitPop Tenn_nomod Tenn_pulsed
do
    for pop in EAS EUR nonAfr AFR modHum
    do
        echo using model $mod with population $pop
        ./get_samples.sh $mod $pop $len &
        cd ../../../old_v/msprime_scripts
        ./get_samples.sh $mod $pop $len &
        cd - > /dev/null
        wait

        ./comp_samples.sh $mod $pop
        echo no diffs...
    done
done
echo done

