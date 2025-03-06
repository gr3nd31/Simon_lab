#! /bin/bash

seq_file=$1

# Creates blast database in the 'reference' directory
echo "Creating blast db..."
cd references

# Converts T's to U's and concats references
for i in $(ls);
do
    cat $i | sed 's/U/T/g' >> db.fasta
done
# Makes a combined database
makeblastdb -in db.fasta -dbtype nucl
cd ../

# Runs a preliminary alignment
echo "Running alignment..."
blastn -db references/db.fasta -query $seq_file -outfmt 6 -out db_aligned.tsv
rm -r references/db.fasta*

# Generates a list of reads that likely belong to each reference
echo "Binning reads..."
Rscript bin/refine.R
rm db_aligned.tsv

echo "Pulling reads from source file..."
echo ""
cd reads
for i in $(ls);
do
    cd $i
    for j in $(ls ../../references/);
        do
            t=$(cat ../../references/$j | head -n 1 | sed 's/>//g')
            if [ ! -f "$t.json" ]; then
                if [ $t == $i ]; then
                    python3 ../../bin/pullbacks.py -l $(ls *txt) -r ../../$1
                    cp ../../references/$j ref.fasta
                    makeblastdb -in ref.fasta -dbtype nucl
                    blastn -query out_reads.fasta -db ref.fasta -outfmt 15 -out $i.json
                    python3 ../../bin/granulate.py -a $i.json -g ref.fasta
                    rm ref.fasta*
                fi

            else
                echo Skipping sample $j
            fi
        done
    cd ../
done
echo ""
echo "Read refinement, alignment, and granulation complete."