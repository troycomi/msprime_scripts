#!/bin/bash
#!/usr/bin/python/2.7.14

#SBATCH --get-user-env
#SBATCH --mem=35G
#SBATCH --time 6-0
#SBATCH --qos=long
#SBATCH --output=slrun.msprime.Sriram.2pulse.0_to_20pc.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

dir=~/SimulatedDemographic/msprime/
echo $dir

for n in 0.0 0.01 0.02 0.03 0.04 0.05 0.07 0.10 0.15 0.2; do
#for n in $(seq 0.00 0.01 0.05); do
	if [ $n = "0.00" ] ; then 
		n=$(echo 0.0)
	elif [ $n = "0.10" ] ; then
		n=$(echo 0.1)
	elif [ $n = "0.20" ] ; then
		n=$(echo 0.2)
	fi

	for d in 0.00 0.01 0.02 0.03 0.04 0.05 0.07 0.10; do
#	for d in $(seq 0.01 0.01 0.01); do
		if [ $d = "0.00" ]; then
			d=$(echo 0.0)
		elif [ $d = "0.10" ] ; then
			d=$(echo 0.1)
		elif [ $d = "0.20" ] ; then
			d=$(echo 0.2)
		fi

		#echo -e Neand_Pulse1: "$n" , Neand_Pulse2: "$d" , Seed: "${SGE_TASK_ID}"
		
#		#echo Tenn
#		time python $dir/bin/msprime.Admixture_Simulate.py \
#			-p nonAfr \
#			-o Tenn \
#			-s ${SLURM_ARRAY_TASK_ID} \
#			-i 2 \
#			-n $n \
#			-d $d \
#			-e 1006 \
#			-a 2040 \
#			-c haplo \
#			-l 1e7
#		
		#echo Sriram
		echo ''
		time python $dir/bin/msprime.Admixture_Simulate.01_14_18.py \
			-p nonAfr \
			-o Sriram \
			-s ${SLURM_ARRAY_TASK_ID} \
			-i 2 \
			-n $n \
			-d $d \
			-e 1006 \
			-a 2040 \
			-c haplo \
			-l 1e7
		
#		#echo SplitPop
#		echo ''
#		time python $dir/bin/msprime.Admixture_Simulate.py \
#			-p nonAfr \
#			-o SplitPop \
#			-s ${SLURM_ARRAY_TASK_ID} \
#			-i 2 \
#			-n $n \
#			-d $d \
#			-e 1006 \
#			-a 2040 \
#			-c haplo \
#			-l 1e7
#
		
	done
done
echo ''
echo FIN
date
