#!/bin/bash

for pop in AFR EUR ASN; do
    awk -v pop=$pop 'NR != 1 { if($2 == pop) print $1}' base.popfile > $pop.pop
    vcftools --gzvcf test.vcf.gz --site-pi --keep $pop.pop --maf 0.05 2> /dev/null
    res=$(awk 'NR != 1 {s += $3} END{print s / NR}' out.sites.pi)
    echo $pop $res
done

vcftools --gzvcf test.vcf.gz --weir-fst-pop AFR.pop --weir-fst-pop EUR.pop --maf 0.05 2> /dev/null
res=$(awk 'NR != 1 {s += $3} END{print s / NR}' out.weir.fst)
echo AFR EUR $res

vcftools --gzvcf test.vcf.gz --weir-fst-pop AFR.pop --weir-fst-pop ASN.pop --maf 0.05 2> /dev/null
res=$(awk 'NR != 1 {s += $3} END{print s / NR}' out.weir.fst)
echo AFR ASN $res

vcftools --gzvcf test.vcf.gz --weir-fst-pop EUR.pop --weir-fst-pop ASN.pop --maf 0.05 2> /dev/null
res=$(awk 'NR != 1 {s += $3} END{print s / NR}' out.weir.fst)
echo EUR ASN $res
