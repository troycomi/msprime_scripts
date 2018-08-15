#$ -S /bin/bash
#$ -l mfree=3G

for a in $(seq 0.0 0.01 0.2); do
	python msprime.Rajiv.AWedit.py -p nonAfr -o corr -s 9 -i 1 -a $a
	
	~/AkeyRotation/bin/AdmixTools/bin/qpF4ratio -p parfile.F4stat.sim_"$a"pc > sim.output_F4stat_"$a"pc
	~/AkeyRotation/bin/AdmixTools/bin/qpDstat -p parfile.Dstat.sim_"$a"pc > sim.output_Dstat_"$a"pc
	
	pct_int=$(cat corrnonAfr[0-9]_"$a"/*.bed | awk 'BEGIN {OFS="\t"} {len_bp=$3-$2 ; sum_len_bp+=len_bp} END {print sum_len_bp/20/2000/1000000}')
	echo -e $a'\t'$pct_int >> corr_pct_int.txt
done
