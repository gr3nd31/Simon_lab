# How to Use

## These scripts use the ViennaRNA package to:

+ Randomly generate new hairpins with certain parameters, and record the thermodynamic properties of these hairpins (**make_a_hairpin.py**). The following parameters can be called:
    -   -n, --NUMBEROFITERATIONS: The number of hairpin variants to generate. (Integer: Default is 1)
    -   -s, --SEQUENCE: Path to a text file containing an RNA sequence to be placed on the 5' side of the hairpin (String).
    -   -S, --STARTSEQUENCE: The relative position where a given sequence starts on the 5' side of the hairpin. (Float: 0.0 - 1.0. The default is 0, meaning the sequence is placed at the base of the hairpin. A position of 0.5 will place the sequence in the middle of the hairpin, while a position of 1 will place the sequence near the apical loop).
    -   -b, --BULGETYPE: Type of bulge to place in the hairpin (String: `symmetrical`, `left`, or `right`). By default, no bulges are placed in the hairpin.
    -   -B, --BULGEPOSITION: The relative position where bulge is placed. (Float: 0.0 - 1.0. The default is 0.5, meaning the bulge is placed in the middle of the hairpin, while a position of 1 will place the bulge near the apical loop and a position of 0 places the bulge at the base of the hairpin).
    -   -c, --BULGESIZE: Number of nucleotides within the bulge (Integer: Default is 5).
    -   -R, --REPEATS: Flags repeats of *n* length, where *n* is the given integer (Integer: No default).
    -   -D, --DETAILS: Report the apical, bulge, and repeat sequences rather than the number.
    -   -p, --PAIREDPERCENT:Percent of the hairpin that is paired (Float: 0.0 - 1.0. Default is 1.0, meaning the 5' side is 100% complementary to the 3' side).
    -   -l, --LENGTH: Length of the stem (Integer: Default is 30 bases).
    -   -a, --APICALLOOP: Size of the apical loop (Integer: Defaultis 5 bases).
    -   -A, --APICALSEQUENCE: Specific sequence to be placed in the apical loop (String: No default).
    -   -o, --OUT: Name of the output CSV file. If file already exists, the data is appended to the existing CSV file (String: Default is 'data.csv').

+ Take a FASTA file of sequences and return the ViennaRNA folded structures and thermodynamic properties (**give_a_hairpin.py**: `python give_a_hairpin.py -s SEQUENCES.fasta`). The following parameters can be called:
    -   -s, --SEQUENCE: Path to a fasta file containing an RNA sequence(s) to be folded and parsed (String).
    -   -l, --LENGTH: Maximum length to be folded. If a given RNA sequence is longer than this parameter, it will be split into smaller segments for folding.
    -   -R, --REPEATS: Flags repeats of *n* length, where *n* is the given integer (Integer: No default).
    -   -D, --DETAILS: Report the apical, bulge, and repeat sequences rather than the number.
    -   -o, --OUT: Name of the output CSV file. If file already exists, the data is appended to the existing CSV file (String: Default is 'data.csv').

+ Search a larger dotbracket structure and sequence for a substructure and return the RNA sequence(s) that correspond to that sequence in a FASTA format (**find_a_hairpon.py**: `python find_a_hairpin.py -s SEQUENCE.txt`). The following parameters can be called:
    -   -s, --SEQUENCE: Path to a fasta file containing an RNA sequence and dot-bracket structure. In this file, the first line should be the seuqence name, the second line is RNA sequence, and the third line is the dot-bracket structure (String).
    -   -m, --MINLENGTH: The minimum length a structure must be in order to be recorded (Integer: Default is 11 bases).
    -   -f, --FIND:  Regular expression for the structure to be found using 'C' for '(', 'G' for ')', and 'U' for '.' . Default is a simple hairpin 'C+C+[CU]+C+U+G+[GU]+G+G+'. A YSS would be 'C+[CU]+U+G+[GU]+C+[CU]+G+[GU]+G+'.
    -   -p, --PREFIX: Prefix to be added to the subsequence name. (String: Default name is derived from the first line of the text file).
    -   -o, --OUT: Name of the output fasta file (String: Default is 'sequences.fasta').

+ Take a given sequence and mutate the sequence to generate alternative RNA sequences (and maintain coding) with differing thermodynamic properties (**mutate_a_hairpon.py**: `python mutate_a_hairpin.py -s ORIGINAHAIRPIN.txt -c STARTING_BASE_NUMBER -n NUMBER_OF_ITERATIONS`). The following parameters can be called:
    -   -s, --SEQUENCE: Path to a fasta file containing a sequence to be mutated (String).
    -   -p, --PERCENTDIFFERENT:Percent of the sequence that is should be changed (Float: 0.0 - 1.0. Default is 0.5, meaning that 50% of the possible changes will be made).
    -   -R, --RANDOMIZEPERCENT: If given, the PERCENTDIFFERENT will be randomly selected from 1-100% for each iteration.
    -   -c, --CODING: If given, this parameter (Integer) will change codons to keep amino acid coding rather than changing individual bases. The number should correspond to the sequence position where the coding begins.
    -   -e, --ENDCODING: Position (Integer) at which coding should be ended if the -c parameter has been given.
    -   -u, --USAGETABLE: Path to the codon usage table. If no table is given, codon changes are random. If a table is given, codon selection is weighted based on the given table (String: Several tables are present in the bin/codon_tables directory).
    -   -n, --NUMBEROFITERATIONS: The number of variants to generate (Integer: Default is 100).
    -   o, --OUT: Name of the output fasta file (String).

Additionally, there are additional scripts to assist in analyses:

+ **LbG_graph.R** takes a output csv file from **give_a_hairpin.py**/**give_a_hairpin.py** and graphs a Length by MFE dotplot containing the linear equation as a title. (Run from command line: `Rscript LbG_graph.R data.csv`)
+ **pe_slope.R** takes a dataframe as an argument and calculates the slope of Positional Entropy across the given structure. (Run within R: `datum <- pe_slope(datum)`)
+ **check_codons.py** takes a coding sequence of RNA can checks the codon usage against a given codon table. The ouput is a dot plot of expected frequency vs the found codon usage.

Please see code (or JMN) for specific usage.
