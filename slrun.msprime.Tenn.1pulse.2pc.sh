#!/bin/bash
#!/usr/bin/python/2.7.14

#SBATCH --get-user-env
#SBATCH --mem=75G
#SBATCH --time 100-0
#SBATCH --qos=long
#SBATCH --output=slrun.msprime.Tenn.1pulse.2pc.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

dir=~/SimulatedDemographic/msprime/
echo $dir

time python $dir/bin/msprime.Admixture_Simulate.py \
	-p nonAfr \
	-o Tenn \
	-s ${SLURM_ARRAY_TASK_ID} \
	-i 1 \
	-n 0.02 \
	-c F4Dstat \
	-e 500000 \
	-a 1 \
	-l 100e6

echo ''
echo FIN
date
