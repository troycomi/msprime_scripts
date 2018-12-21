#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mem=8G
#SBATCH --time 4-0
#SBATCH --qos=long
#SBATCH --output=slrun.merge_bed.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
dir=/Genomics/akeylab/abwolf/SimulatedDemographic/msprime/
mdl=$(echo Tenn)
echo $mdl
echo ''

for rep in $(seq 48 1 50); do
	for n1 in $(seq 0.0 0.01 0.2); do
		if [ $n1 = "0.00" ] ; then 
			n1=$(echo 0.0)
		elif [ $n1 = "0.10" ] ; then
			n1=$(echo 0.1)
		elif [ $n1 = "0.20" ] ; then
			n1=$(echo 0.2)
		fi
		
		for n2 in $(seq 0.0 0.01 0.1); do
			if [ $n2 = "0.00" ] ; then 
				n2=$(echo 0.0)
			elif [ $n2 = "0.10" ] ; then
				n2=$(echo 0.1)
			elif [ $n2 = "0.20" ] ; then
				n2=$(echo 0.2)
			fi
			
			simname=$( echo "$mdl"_nonAfr_"$rep"_n1_"$n1"_n2_"$n2".bed )
			#echo $simname
			
			if [ -e $dir/$mdl/$simname.gz ] && [ ! -e $dir/$mdl/$simname.merged.gz ]; then
				echo $dir/$mdl/$simname.gz
				zcat $dir/$mdl/$simname.gz \
					| awk 'BEGIN {OFS="\t"} {if($2<$3) print $0}' \
					| sort-bed - \
					| bedops --merge - \
					| awk 'BEGIN {OFS="\t"} {print "'$rep'", $2, $3, $1}' \
					| gzip -c - \
					> $dir/$mdl/$simname.merged.gz
			else
				echo MISSING: $simname
			fi

		done
	done
done

echo ''
echo fin
date

