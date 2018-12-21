#!/bin/bash

oth=../../../old_v/msprime_scripts/output
here=version_1

set -e
mod="$1"
pop="$2"

cmp $here/debug.txt $oth/debug.out
cmp $here/${mod}_${pop}_2_n1_0.1_n2_0.1.popfile $oth/${mod}.popfile

#skip version difference on line 2
cmp -i 43 <(zcat $here/${mod}_${pop}_2_n1_0.1_n2_0.1.vcf.gz) $oth/vcf.out
cmp <(zcat $here/${mod}_${pop}_2_n1_0.1_n2_0.1.bed.merged.gz) $oth/haplo.out
cmp <(zcat $here/${mod}_${pop}_2_n1_0.1_n2_0.1.ils.bed.merged.gz) $oth/ils.out
cmp <(zcat $here/${mod}.eigenstratgeno.n1_0.1_n2_0.1_t_145_2.gz) <(zcat $oth/${mod}.eigenstratgeno.n1_0.1_n2_0.1_t_145_mAfB_0.00014_mBAf_0.00015_mAfEu_1.5e-05_mEuAf_2.5e-05_2.gz)
cmp <(zcat $here/${mod}.snp.n1_0.1_n2_0.1_t_145_2.gz) <(zcat $oth/${mod}.snp.n1_0.1_n2_0.1_t_145_mAfB_0.00014_mBAf_0.00015_mAfEu_1.5e-05_mEuAf_2.5e-05_2.gz)
cmp <(zcat $here/${mod}.ind.n1_0.1_n2_0.1_t_145_2.gz) <(zcat $oth/${mod}.ind.n1_0.1_n2_0.1_t_145_mAfB_0.00014_mBAf_0.00015_mAfEu_1.5e-05_mEuAf_2.5e-05_2.gz)
