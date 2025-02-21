#! /bin/bash

seq_file=$1

# Creates blast database in the 'reference' directory
echo "Creating blast db..."
cd references

for i in $(ls);
do
    cat $i | sed 's/U/T/g' >> db.fasta
done
makeblastdb -in db.fasta -dbtype nucl
cd ../

echo "Running alignment..."
blastn -db references/db.fasta -query $seq_file -outfmt 6 -out db_aligned.tsv
rm -r references/db.fasta*

echo "Binning reads..."
Rscript bin/refine.R
rm db_aligned.tsv

echo "Pulling reads from source file..."
cd reads
for i in $(ls);
do
    cd $i
    #echo $i
    python3 ../../bin/pullbacks.py -l $(ls *txt) -r ../../$1
    for j in $(ls ../../references/);
    do
         t=$(cat ../../references/$j | head -n 1 | sed 's/>//g')
         #echo $t
         if [ $t == $i ]; then
            cp ../../references/$j ref.fasta
            makeblastdb -in ref.fasta -dbtype nucl
            blastn -query out_reads.fasta -db ref.fasta -outfmt 15 -out $i.json
            python3 ../../bin/granulate.py -a $i.json -g ref.fasta
            rm ref.fasta*
        fi
    done
    cd ../
done

echo "Read refinement, alignment, and granulation complete."