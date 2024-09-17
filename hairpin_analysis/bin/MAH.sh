#! /bin/bash

mkdir bulge_data/$1
for i in 2 3 4 5 6 ; 
    do
        python3 bin/make_a_hairpin.py -o bulge_data/$1/data.csv -n 100 -b symmetric -c $i -l $1 -B 0.2 -p 0.9
    done;