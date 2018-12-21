#$ -S /bin/bash
#$ -l mfree=3G

for n in $(seq 0.0 0.01 0.2); do
	if [ $n = "0.00" ] ; then 
		n=$(echo 0.0)
	elif [ $n = "0.10" ] ; then
		n=$(echo 0.1)
	elif [ $n = "0.20" ] ; then
		n=$(echo 0.2)
	fi
	for d in $(seq 0.0 0.01 0.2); do
		if [ $d = "0.00" ]; then
			d=$(echo 0.0)
		elif [ $d = "0.10" ] ; then
			d=$(echo 0.1)
		elif [ $d = "0.20" ] ; then
			d=$(echo 0.2)
		fi
		echo -e Neand_Pulse1: "$n" , Neand_Pulse2: "$d" , Seed: "${SGE_TASK_ID}"
	#	python msprime.Tennessen.py -p nonAfr -o Tenn -s ${SGE_TASK_ID} -i 2 -n $n -d $d
	#	python msprime.Sriram.py -p nonAfr -o Sriram -s ${SGE_TASK_ID} -i 2 -n $n -d $d
	#	python msprime.SplitPop.py -p nonAfr -o Split -s ${SGE_TASK_ID} -i 2 -n $n -d $d
		
		~/AkeyRotation/bin/AdmixTools/bin/qpF4ratio -p parfile.F4stat.Split.n1_"$n"_n2_"$d"_1 > Split.output_F4stat.n1_"$n"_n2_"$d"_1
		~/AkeyRotation/bin/AdmixTools/bin/qpF4ratio -p parfile.F4stat.Tenn.n1_"$n"_n2_"$d"_1 > Tenn.output_F4stat.n1_"$n"_n2_"$d"_1
		~/AkeyRotation/bin/AdmixTools/bin/qpF4ratio -p parfile.F4stat.Sriram.n1_"$n"_n2_"$d"_1 > Sriram.output_F4stat.n1_"$n"_n2_"$d"_1
	done
done

echo fin
