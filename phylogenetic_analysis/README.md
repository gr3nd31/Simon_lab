# Dreamcatcher - short script for generating phylogenetic trees 

This bash script takes a text file of Genbank accession numbers and returns a raxml tree. First it will check if the reference sequence is already available and, if not, it will download the file locally and append the number and name to a database. Furthermore, if the accession number is followed by `,X,Y` where X and Y are integers, it will pull only the sequence from X to Y for alignment. This allows you to quickly run phylogeny on designated sequences. Alternatively, if the nucleotide database contains ORF annotation (Format: ORF_NAME1[START:STOP], ORF_NAME2[START:STOP]...), then an orf name can be called instead. This is currently set up to run only nucleotides but may be extended to amino acids in the future.

## Usage

bash PATH/TO/dreamcatcher.sh full/fast nucl LIST.TXT ORF_NAME
