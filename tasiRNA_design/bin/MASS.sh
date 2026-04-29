#! /bin/bash

#Step 1: pull the siRNA sequences from the siRNA file
python3 splitIDs.py -s $1

echo "Running siRNA alignment..."
makeblastdb -in $2 -dbtype nucl -parse_seqids
blastn -db $2 -query out.fasta -outfmt 6 -out db_aligned.tsv -task blastn
rm out.fasta
rm $2.n*
echo "siRNA alignment complete. Gathering data..."
#python3 MASS_ID.py -s $1 -g $2
Rscript MASSive.R db_aligned.tsv $2 $1
echo "siRNA with at least 50% coverage in given sequences identified."
echo "Complete."