#!/bin/bash
#!/usr/bin/python/2.7.14

#SBATCH --get-user-env
#SBATCH --mem=1G
#SBATCH --time 10:00:00
#SBATCH --qos=1day
#SBATCH --output=slrun.F4stat.concat_alpha.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

mdl=$1
dir=~/SimulatedDemographic/msprime/

echo $mdl
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
	
		for i in $(seq 1 1 100); do
			if [ -e F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_"$i".gz ]; then
				zcat F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_"$i".gz \
					| tail -n 5 \
					| head -n 4 \
					| awk 'BEGIN {OFS="\t"} {print $0 }' \
					>> F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_ALL
			
			elif [ ! -e F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_"$i".gz ]; then
				echo ''
				echo WARNING: F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_"$i".gz does not exist
			fi		
		done
		
		time cat F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_ALL \
			| awk 'BEGIN {OFS="\t"} {print $0, "'$n'", "'$d'"}' \
			| gzip -c > F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_ALL.gz
		rm F4stat."$mdl".n1_"$n"_n2_"$d"_t_350_ALL
	done
done
echo ''
echo FIN
date
