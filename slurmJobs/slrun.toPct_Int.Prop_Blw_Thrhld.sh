#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mem=25G
#SBATCH --time 6-0
#SBATCH --qos=1wk
#SBATCH --output=slrun.toPct_Int.Prop_Blw_Thrhld.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
dir=/Genomics/akeylab/abwolf/SimulatedDemographic/msprime/
mdl=$( echo $1 )
echo $mdl
echo ''

pop_size=3046
thrhld=$(awk 'BEGIN {print 0.000316}')
number_chrom=500


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
	
		file=$(echo "$mdl"_nonAfr_ALL_n1_"$n1"_n2_"$n2".bed.merged.5_to_10Mb)
		if [ -e $dir/$mdl/$file.gz ]; then
			echo ''
			echo n1=$n1 n2=$n2    $file.gz
			echo ''
			zcat $dir/$mdl/$file.gz \
				| awk 'BEGIN {OFS="\t"} {print $0,"'$n1'","'$n2'"}' \
				> $dir/$mdl/$file

			Rscript ~/SimulatedDemographic/msprime/bin/toPct_Int.1_to_15Mb.R \
				$dir/$mdl/$file \
				$pop_size \
				$mdl \
				$thrhld \
				$number_chrom \
				> $dir/$mdl/$file.prop_blw_thrhld
			rm $dir/$mdl/$file
		elif [ ! -e $dir/$mdl/$file.gz ]; then
			echo ''
			echo WARNING: $file.gz does not exist
			echo''
		fi
	done
done
echo ''
echo FIN
date
