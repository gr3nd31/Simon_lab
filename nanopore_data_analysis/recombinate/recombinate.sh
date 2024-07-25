#! /bin/bash

db_file=TCV_et_al.fasta
read_file=tcv_reads.fasta
window_length=80
window_overlap=0
prefix=""
blast_task="blastn"

#-------------------------

# Creates blast database in the 'reference' directory
echo "Creating blast db..."
if [ ! -d reference ]; then
    mkdir reference
fi
cp $db_file reference
makeblastdb -in reference/$db_file -dbtype nucl

# Split the 
echo "Splitting reads into windows..."
python windowfy.py -r $read_file -w $window_length

echo "Beginning blast alignments..."
for dir in recombo/*/
do
    target=${dir}windowed.fasta
    dips=${dir}windowed.json
    blastn -db reference/$db_file -query $target -outfmt 15 -out $dips -task $blast_task
    python recombo.py -i $dips
    rm ${dir}windowed.json
done
echo "Recombination check complete. Now analyze with R script"