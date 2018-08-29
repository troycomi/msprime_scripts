#!/bin/bash
#!/usr/bin/python/2.7.14

#SBATCH --get-user-env
#SBATCH --mem=2G
#SBATCH --time 00:30:00
#SBATCH --qos=1hr
#SBATCH --output=slrun.F4stat.Tenn_Sriram_Split.1pulse.migration.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

dir=~/SimulatedDemographic/msprime/
echo $dir

n=$( echo 0.0 )
d=$( echo 0.0 )
eur=$( echo 1006 )
asn=$( echo 1008 )
afr=$( echo 216 )
len=$( echo 1e5 )

m_AF_B=$( echo 0.0 )
m_B_AF=$( echo 0.0 )
m_AF_EU=$( echo 0.0 )
m_EU_AF=$( echo 0.1 )




time python $dir/bin/msprime.Admixture_Simulate.py \
	-p nonAfr \
	-o Tenn \
	-s ${SLURM_ARRAY_TASK_ID} \
	-i 2 \
	-n $n \
	-d $d \
	-e $eur \
	-a $asn \
	-r $afr \
	--migration_AF_B $m_AF_B \
	--migration_B_AF $m_B_AF \
	--migration_AF_EU $m_AF_EU \
	--migration_EU_AF $m_EU_AF \
	--migration_AF_AS 0.0 \
	--migration_AS_AF 0.0 \
	--migration_EU_AS 0.0 \
	--migration_AS_EU 0.0 \
	-c F4Dstat \
	-l $len

time gunzip \
	parfile.F4stat.Tenn.n1_"$n"_n2_"$d"_t_350_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF"_${SLURM_ARRAY_TASK_ID}.gz \
	Tenn.eigenstratgeno.n1_"$n"_n2_"$d"_t_350_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF"_${SLURM_ARRAY_TASK_ID}.gz \
	Tenn.snp.n1_"$n"_n2_"$d"_t_350_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF"_${SLURM_ARRAY_TASK_ID}.gz \
	Tenn.ind.n1_"$n"_n2_"$d"_t_350_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF"_${SLURM_ARRAY_TASK_ID}.gz



time qpF4ratio -p parfile.F4stat.Tenn.n1_"$n"_n2_"$d"_t_350_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF"_${SLURM_ARRAY_TASK_ID} \
	| gzip -c - \
	> F4stat.Tenn.n1_"$n"_n2_"$d"_t_350_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF"_${SLURM_ARRAY_TASK_ID}.gz

# rm parfile.F4stat.Tenn.n1_"$n"_n2_"$d"_t_350_${SLURM_ARRAY_TASK_ID} \
# 	Tenn.eigenstratgeno.n1_"$n"_n2_"$d"_t_350_${SLURM_ARRAY_TASK_ID} \
# 	Tenn.snp.n1_"$n"_n2_"$d"_t_350_${SLURM_ARRAY_TASK_ID} \
# 	Tenn.ind.n1_"$n"_n2_"$d"_t_350_${SLURM_ARRAY_TASK_ID}


#		#echo Sriram
#		echo ''
#		time python $dir/bin/msprime.Admixture_Simulate.py \
#			-p nonAfr \
#			-o Sriram \
#			-s ${SLURM_ARRAY_TASK_ID} \
#			-i 2 \
#			-n $n \
#			-d $d \
#			-c haplo \
#			-l 1e7

#		#echo SplitPop
#		echo ''
#		time python $dir/bin/msprime.Admixture_Simulate.py \
#			-p nonAfr \
#			-o SplitPop \
#			-s ${SLURM_ARRAY_TASK_ID} \
#			-i 2 \
#			-n $n \
#			-d $d \
#			-c haplo \
#			-l 1e7
#
echo ''
echo FIN
date
