#!/bin/bash

diff Tenn.popfile v1/Tenn.popfile | head

for f in Tenn*.gz; do
    echo $f
    diff <(zcat $f) <(zcat v1/$f) | head
done

f=parfile.F4stat.Tenn.n1_0.1_n2_0.1_t_350_2.gz
echo $f
diff <(zcat $f) <(zcat v1/$f) | head

f=Tenn.popfile
echo $f
