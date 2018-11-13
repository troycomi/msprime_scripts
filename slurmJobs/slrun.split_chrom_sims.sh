#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mem=25G
#SBATCH --time 10-0
#SBATCH --qos=long
#SBATCH --output=slrun.split_chrom_sims.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
dir=/Genomics/akeylab/abwolf/SimulatedDemographic/msprime/
mdl=$( echo $1 )
echo $mdl
echo ''

for n1 in 0.0 0.01 0.02 0.03 0.04 0.05 0.07 0.10 0.15 0.2; do
#for n1 in $(seq 0.0 0.01 0.2); do
	if [ $n1 = "0.00" ] ; then 
		n1=$(echo 0.0)
	elif [ $n1 = "0.10" ] ; then
		n1=$(echo 0.1)
	elif [ $n1 = "0.20" ] ; then
		n1=$(echo 0.2)
	fi
	
	for n2 in 0.00 0.01 0.02 0.03 0.04 0.05 0.07 0.10; do
	#for n2 in $(seq 0.0 0.01 0.1); do
		if [ $n2 = "0.00" ] ; then 
			n2=$(echo 0.0)
		elif [ $n2 = "0.10" ] ; then
			n2=$(echo 0.1)
		elif [ $n2 = "0.20" ] ; then
			n2=$(echo 0.2)
		fi
		
		if [ -e $dir/$mdl/"$mdl"_nonAfr_1_n1_"$n1"_n2_"$n2".bed.merged.gz ]; then
			echo $dir/$mdl/"$mdl"_nonAfr_1_n1_"$n1"_n2_"$n2".bed.merged.gz
		
			time cat \
				$dir/$mdl/"$mdl"_nonAfr_*_n1_"$n1"_n2_"$n2".bed.merged.gz \
				> $dir/$mdl/"$mdl"_nonAfr_ALL_n1_"$n1"_n2_"$n2".bed.merged.gz 
			
			time python \
				~/SimulatedDemographic/msprime//bin/split_chromosomes.temp.loop.py \
				$dir/$mdl/"$mdl"_nonAfr_ALL_n1_"$n1"_n2_"$n2".bed.merged.gz \
				$dir/$mdl/"$mdl"_nonAfr_ALL_n1_"$n1"_n2_"$n2".bed.merged.5_to_10Mb.gz \
				10000000
			echo ''
		elif [ ! -e $dir/$mdl/"$mdl"_nonAfr_1_n1_"$n1"_n2_"$n2".bed.merged.gz ]; then
			echo WARNING: MISSING 	n1=$n1 n2=$n2
			echo ''
		fi
	done
done

echo ''
echo FIN
date
