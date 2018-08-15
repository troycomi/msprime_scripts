#$ -S /bin/bash
#$ -l mfree=3G

echo -e Model'\t'Admix_Prop'\t'pct_int'\t'F_EUR_N1'\t'FZ_EUR_N1'\t'F_ASN_N1'\t'FZ_ASN_N1'\t'F_EUR_N2'\t'FZ_EUR_N2'\t'F_ASN_N2'\t'FZ_ASN_N2'\t'D_EUR_N1'\t'DZ_EUR_N1'\t'D_ASN_N1'\t'DZ_ASN_N1'\t'D_EUR_N2'\t'DZ_EUR_N2'\t'D_ASN_N2'\t'DZ_ASN_N2 >> PctInt_F4_D.txt

for sim in Tenn Sriram Split; do
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
		
			pct_int=$(echo NA)

			F4folder=$(echo "$sim".output_F4stat)	
			F4file=$(echo "$sim".output_F4stat.n1_"$n"_n2_"$d"_1)
			Dfolder=$(echo "$sim".output_Dstat)
			Dfile=$(echo "$sim".output_Dstat.n1_"$n"_n2_"$d"_1)

			F_EUR_N1=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $11}' )
			F_EUR_N2=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $11}' )
			F_ASN_N1=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $11}' )
			F_ASN_N2=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $11}' )
			
			FZ_EUR_N1=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $13}' )
			FZ_EUR_N2=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $13}' )
			FZ_ASN_N1=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $13}' )
			FZ_ASN_N2=$( cat $F4folder/$F4file | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $13}' )
		
		
			D_EUR_N1=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $6}' )
			D_EUR_N2=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $6}' )
			D_ASN_N1=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $6}' )
			D_ASN_N2=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $6}' )
		
			DZ_EUR_N1=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $7}' )
			DZ_EUR_N2=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $7}' )
			DZ_ASN_N1=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $7}' )
			DZ_ASN_N2=$( cat $Dfolder/$Dfile | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $7}' )
		
			echo -e $sim'\t'n1_"$n"_n2_"$d"'\t'$pct_int'\t'$F_EUR_N1'\t'$FZ_EUR_N1'\t'$F_ASN_N1'\t'$FZ_ASN_N1'\t'$F_EUR_N2'\t'$FZ_EUR_N2'\t'$F_ASN_N2'\t'$FZ_ASN_N2'\t'$D_EUR_N1'\t'$DZ_EUR_N1'\t'$D_ASN_N1'\t'$DZ_ASN_N1'\t'$D_EUR_N2'\t'$DZ_EUR_N2'\t'$D_ASN_N2'\t'$DZ_ASN_N2 >> PctInt_F4_D.txt
		done
	done
done
echo FIN


