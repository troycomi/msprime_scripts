#!/bin/bash
source ~/.bashrc

s=${RANDOM}
#s=1000

cmd=$( echo " python ~/SimulatedDemographic/msprime/bin/msprime.Admixture_Simulate.py \n
        -o Tenn_pulsed \n
        -s $s -i 2 \n
        -n 0.1 -d 0.0 -e 100 -a 100 -r 100 -l 1000000 \n
        --migration_AF_B 0.0 --migration_B_AF 0.0 \n
        --migration_AF_EU 0.0 --migration_EU_AF 0.1 \n
        --migration_AF_AS 0.0 --migration_AS_AF 0.0 \n
        --migration_EU_AS 0.0 --migration_AS_EU 0.0 \n")
echo -e $cmd

for p in modHum nonAfr EUR EAS AFR; do
	echo $p
	cmd2=$( echo -e " $cmd \n
	-p $p \n
	-c haplo | cut -f 4 | sort -n - | uniq - | wc -l ")
	#echo -e $cmd2
	eval $( echo -e $cmd2 )
done
echo ''




#        --migration_AF_B 0.0 --migration_B_AF 0.0 \n
#        --migration_AF_EU 0.0 --migration_EU_AF 0.0 \n
#        --migration_AF_AS 0.0 --migration_AS_AF 0.0 \n
#        --migration_EU_AS 0.0 --migration_AS_EU 0.0 \n
# $(awk 'BEGIN {print 0.1/3000}')
# $(awk 'BEGIN {print 0.1/720}')
