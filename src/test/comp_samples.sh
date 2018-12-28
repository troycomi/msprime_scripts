#!/bin/bash

set -e
mod="$1"
pop="$2"

oth=version_1/${3}_${pop}_${mod}
here=version_2/${3}_${pop}_${mod}

cmp $here/debug.txt $oth/debug.txt
cmp $here/${mod}_${pop}_2_n1_0.1_n2_0.1.popfile $oth/${mod}_${pop}_2_n1_0.1_n2_0.1.popfile

cmp <(zcat $here/${mod}_${pop}_2_n1_0.1_n2_0.1.vcf.gz) <(zcat $oth/${mod}_${pop}_2_n1_0.1_n2_0.1.vcf.gz)
cmp <(zcat $here/${mod}_${pop}_2_n1_0.1_n2_0.1.bed.merged.gz) <(zcat $oth/${mod}_${pop}_2_n1_0.1_n2_0.1.bed.merged.gz)
cmp <(zcat $here/${mod}_${pop}_2_n1_0.1_n2_0.1.ils.bed.merged.gz) <(zcat $oth/${mod}_${pop}_2_n1_0.1_n2_0.1.ils.bed.merged.gz)
cmp <(zcat $here/${mod}.eigenstratgeno.n1_0.1_n2_0.1_t_145_2.gz) <(zcat $oth/${mod}.eigenstratgeno.n1_0.1_n2_0.1_t_145_2.gz)
cmp <(zcat $here/${mod}.snp.n1_0.1_n2_0.1_t_145_2.gz) <(zcat $oth/${mod}.snp.n1_0.1_n2_0.1_t_145_2.gz) 
cmp <(zcat $here/${mod}.ind.n1_0.1_n2_0.1_t_145_2.gz) <(zcat $oth/${mod}.ind.n1_0.1_n2_0.1_t_145_2.gz) 
