#$ -S /bin/bash
#$ -l mfree=5G


for n in $(seq 0.00 0.01 0.20); do
        if [ $n = "0.00" ] ; then
               	n=$(echo 0.0)
        elif [ $n = "0.10" ] ; then
               	n=$(echo 0.1)
        elif [ $n = "0.20" ] ; then
               	n=$(echo 0.2)
        fi

        for d in $(seq 0.00 0.01 0.20); do
                if [ $d = "0.00" ]; then
                       	d=$(echo 0.0)
                elif [ $d = "0.10" ] ; then
                       	d=$(echo 0.1)
                elif [ $d = "0.20" ] ; then
                       	d=$(echo 0.2)
                fi
			
		for s in $(seq 1 1 10); do
			
			echo -e Neand_Pulse1: "$n" , Neand_Pulse2: "$d" , Seed: "$s"

			cat Tenn_nonAfr_"$s"_n1_"$n"_n2_"$d"/Tenn_nonAfr_"$s"_n1_"$n"_n2_"$d".bed \
			| awk 'BEGIN {OFS="\t"} {print '"$s"', $2, $3, $4, '"$n"','"$d"'}' \
			>> Tenn_nonAfr_ALL_n1_"$n"_n2_"$d".bed

                        cat Sriram_nonAfr_"$s"_n1_"$n"_n2_"$d"/Sriram_nonAfr_"$s"_n1_"$n"_n2_"$d".bed \
                        | awk 'BEGIN {OFS="\t"} {print '"$s"', $2, $3, $4, '"$n"','"$d"'}' \
                        >> Sriram_nonAfr_ALL_n1_"$n"_n2_"$d".bed

                        cat Split_nonAfr_"$s"_n1_"$n"_n2_"$d"/Split_nonAfr_"$s"_n1_"$n"_n2_"$d".bed \
                        | awk 'BEGIN {OFS="\t"} {print '"$s"', $2, $3, $4, '"$n"','"$d"'}' \
                        >> Split_nonAfr_ALL_n1_"$n"_n2_"$d".bed
		done
	done
done

echo FIN
