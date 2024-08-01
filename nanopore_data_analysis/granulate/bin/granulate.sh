#! /bin/bash

seq_file=$1
db_file=$2
prefix=$3

# Creates blast database in the 'reference' directory
echo "Creating blast db..."
if [ ! -d reference ]; then
    mkdir reference
fi
cp $db_file reference
makeblastdb -in reference/$db_file -dbtype nucl

echo "Running alignment..."
blastn -db reference/$db_file -query $seq_file -outfmt 15 -out "$prefix""blast.json"
python3 bin/granulate.py -g reference/$db_file -a "$3""blast.json" -t $prefix

rm -r reference