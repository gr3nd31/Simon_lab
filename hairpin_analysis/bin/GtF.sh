#! /bin/bash

cd subs
seqs=$(ls *seq)
pwd
for seq in $seqs; do
    echo $seq
    x=$(echo $seq | sed -r 's/seq/str/g')
    y=$(echo $seq | sed -r 's/seq//g')
    python3 ../bin/find_a_hairpin.py -s $seq -d $x -m 11 -o ../hairpins.fasta -p $y
done

cd ../