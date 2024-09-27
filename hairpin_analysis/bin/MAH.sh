#! /bin/bash

#mkdir bulge_data/$1
for i in 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 ; 
    do
        echo "Generating hairpins with complimentarity of $i"
        for j in 15 30 45 60 75 90
            echo "Generating hairpins with length of $j"
            do
                python3 bin/make_a_hairpin.py -o data_none.csv -n 400 -l $j -p $i
            done;
    done;