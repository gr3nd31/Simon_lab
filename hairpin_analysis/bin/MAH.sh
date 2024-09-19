#! /bin/bash

#mkdir bulge_data/$1
for i in 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 ; 
    do
        python3 bin/make_a_hairpin.py -o data_15_left.csv -n 400 -l 15 -p $i -b left
    done;