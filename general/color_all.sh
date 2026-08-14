#! /bin/bash

for fileName in $(ls $1); do
            bim=$(echo $fileName | sed 's/txt/rnacanvas/')
            echo $fileName
            echo $bim
            python3 ~/Documents/Github/Simon_lab/general/RNAConvictIt.py -p $1/$fileName -r $2/$bim
        done