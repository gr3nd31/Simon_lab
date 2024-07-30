# Dreamcatcher - short script for generating phylogenetic trees
## This bash script takes a text file of Genbank accession numbers and returns a raxml tree. First it will check if the reference sequence is already available and, if not, it will download the file locally and append the number and name to a database. Furthermore, if the accession number is followed by `,X,Y` where X and Y are integers, it will pull only the sequence from X to Y for alignment. This allows you to quickly run phylogeny on designated sequences. Eventually the script will pull these positions automatically from the ORF data in the database but for now these must be entered in manually:

## Usage

bash PATH/TO/dreamcatcher.sh full/fast nucl/aa LIST.TXT
