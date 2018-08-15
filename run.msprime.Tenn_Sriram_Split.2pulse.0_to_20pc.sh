#$ -S /bin/bash
#$ -l mfree=5G

date
echo ''

dir=$(echo ~/AkeyRotation/SimulatedDemographic/msprime/bin/ )

for n in $(seq 0.00 0.01 0.2); do
	
	if [ $n = "0.00" ] ; then 
		n=$(echo 0.0)
	elif [ $n = "0.10" ] ; then
		n=$(echo 0.1)
	elif [ $n = "0.20" ] ; then
		n=$(echo 0.2)
	fi
	
	for d in $(seq 0.00 0.01 0.1); do
		if [ $d = "0.00" ]; then
			d=$(echo 0.0)
		elif [ $d = "0.10" ] ; then
			d=$(echo 0.1)
		elif [ $d = "0.20" ] ; then
			d=$(echo 0.2)
		fi

		#echo -e Neand_Pulse1: "$n" , Neand_Pulse2: "$d" , Seed: "${SGE_TASK_ID}"
		
		#echo Tenn
		time python $dir/msprime.Admixture_Simulate.py -p nonAfr -o Tenn -s ${SGE_TASK_ID} -i 2 -n $n -d $d -c haplo -l 1e7
		#if [ -d Tenn_nonAfr_"${SGE_TASK_ID}"_n1_"$n"_n2_"$d" ]; then
		#	find Tenn_nonAfr_"${SGE_TASK_ID}"_n1_"$n"_n2_"$d" -maxdepth 1 -type f -name "*.sorted.merged" -delete
		#fi
		
		#echo Sriram
		echo ''
		time python $dir/msprime.Admixture_Simulate.py -p nonAfr -o Sriram -s ${SGE_TASK_ID} -i 2 -n $n -d $d -c haplo -l 1e7
		#if [ -d Sriram_nonAfr_"${SGE_TASK_ID}"_n1_"$n"_n2_"$d" ]; then
		#	find Sriram_nonAfr_"${SGE_TASK_ID}"_n1_"$n"_n2_"$d" -maxdepth 1 -type f -name "*.sorted.merged" -delete
		#fi
		
		#echo SplitPop
		echo ''
		time python $dir/msprime.Admixture_Simulate.py -p nonAfr -o SplitPop -s ${SGE_TASK_ID} -i 2 -n $n -d $d -c haplo -l 1e7
		#if [ -d Split_nonAfr_"${SGE_TASK_ID}"_n1_"$n"_n2_"$d" ]; then
		#	find Split_nonAfr_"${SGE_TASK_ID}"_n1_"$n"_n2_"$d" -maxdepth 1 -type f -name "*.sorted.merged" -delete
		#fi

		#~/AkeyRotation/bin/AdmixTools/bin/qpF4ratio -p parfile.F4stat.sim_"$a"pc > sim.output_F4stat_"$a"pc
		#~/AkeyRotation/bin/AdmixTools/bin/qpDstat -p parfile.Dstat.sim_"$a"pc > sim.output_Dstat_"$a"pc
		
	done
done
echo ''
date
echo ''
echo FIN
