# How to Use

## These scripts use the ViennaRNA package to:

+ Randomly generate new hairpins with certain parameters, and record the thermodynamic properties of these hairpins (**make_a_hairpin.py**)
+ Take a FASTA file of sequences and return the ViennaRNA folded structures and thermodynamic properties (**give_a_hairpin.py**: `python give_a_hairpin.py -s SEQUENCES.fasta`)
+ Search a larger dotbracket structure and sequence for a substructure and return the RNA sequence(s) that correspond to that sequence in a FASTA format (**find_a_hairpon.py**: `python give_a_hairpin.py -s SEQUENCE.fasta -d DOTBRACKET.txt`)
+ Take a given sequence and mutate the sequence to generate alternative RNA sequences (and maintain coding) with differing thermodynamic properties (**mutate_a_hairpon.py**: `python mutate_a_hairpin.py -s ORIGINAHAIRPIN.txt -c STARTING_BASE_NUMBER -n NUMBER_OF_ITERATIONS`)

Additionally, there are 2 R scripts to assist in analyses:

+ **LbG_graph.R** takes a output csv file from **give_a_hairpin.py**/**give_a_hairpin.py** and graphs a Length by MFE dotplot containing the linear equation as a title. (Run from command line: `Rscript LbG_graph.R data.csv`)
+ **pe_slope.R** takes a dataframe as an argument and calculates the slope of Positional Entropy across the given structure. (Run within R: `datum <- pe_slope(datum)`)

Please see code (or JMN) for specific usage.