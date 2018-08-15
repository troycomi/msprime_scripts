#$ -S /bin/bash
#$ -l mfree=3G

echo -e Admix_Prop'\t'pct_int'\t'F_EUR_N1'\t'FZ_EUR_N1'\t'F_ASN_N1'\t'FZ_ASN_N1'\t'F_EUR_N2'\t'FZ_EUR_N2'\t'D_EUR_N1'\t'DZ_EUR_N1'\t'D_ASN_N1'\t'DZ_ASN_N1'\t'D_EUR_N2'\t'DZ_EUR_N2'\t'D_ASN_N2'\t'DZ_ASN_N2 >> PctInt_F4_D.txt

for a in $(seq 0.00 0.01 0.20); do
	echo $a

	pct_int=$( cat corrnonAfr[0-9]_"$a"/*.bed | awk 'BEGIN {OFS="\t"} {len_bp=$3-$2 ; sum_len_bp+=len_bp} END {print sum_len_bp/10/2000/1000000}' )

	F_EUR_N1=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $11}' )
	F_EUR_N2=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $11}' )
	F_ASN_N1=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $11}' )
	F_ASN_N2=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $11}' )
	
	FZ_EUR_N1=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $13}' )
	FZ_EUR_N2=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $13}' )
	FZ_ASN_N1=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $13}' )
	FZ_ASN_N2=$( cat sim.output_F4stat_"$a"pc | tail -n 33 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $13}' )


	D_EUR_N1=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $6}' )
	D_EUR_N2=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $6}' )
	D_ASN_N1=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $6}' )
	D_ASN_N2=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $6}' )

	DZ_EUR_N1=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==1) print $7}' )
	DZ_EUR_N2=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==3) print $7}' )
	DZ_ASN_N1=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==2) print $7}' )
	DZ_ASN_N2=$( cat sim.output_Dstat_"$a"pc | tail -n 12 | awk 'BEGIN {OFS="\t"} {if(NR==4) print $7}' )

	echo -e $a'\t'$pct_int'\t'$F_EUR_N1'\t'$FZ_EUR_N1'\t'$F_ASN_N1'\t'$FZ_ASN_N1'\t'$F_EUR_N2'\t'$FZ_EUR_N2'\t'$D_EUR_N1'\t'$DZ_EUR_N1'\t'$D_ASN_N1'\t'$DZ_ASN_N1'\t'$D_EUR_N2'\t'$DZ_EUR_N2'\t'$D_ASN_N2'\t'$DZ_ASN_N2 >> PctInt_F4_D.txt

done

echo fin
