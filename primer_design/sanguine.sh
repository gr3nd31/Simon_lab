#! /bin/bash

db_file=$2
seq_file=$1

#-------------------------

# Creates blast database in the 'reference' directory
echo "Creating blast db..."
if [ ! -d reference ]; then
    mkdir reference
fi
cp $db_file reference
makeblastdb -in reference/$db_file -dbtype nucl

mkdir sequencing_results
cp $seq_file sequencing_results
cd sequencing_results
unzip $seq_file
cd ../

sed -i -e '$a\' sequencing_results/*.seq
cat sequencing_results/*.seq > all_seqs.fasta
rm -r sequencing_results

blastn -db reference/$db_file -query all_seqs.fasta -outfmt 10 -out aligned_seqs.csv
rm all_seqs.fasta
rm -r reference
sed -i -e '1i"Read","Alignment","ID_percent","ID","Gap","U1","Q_start","Q_end","H_start","H_end","Expect","Score"' aligned_seqs.csv
echo "Alignment completed"