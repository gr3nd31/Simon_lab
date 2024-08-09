# Using dreamcatch

The purpose of this script is to take a list of NCBI accession numbers, pull the files, run alignments, and generate a phylogenetic tree. The sequences are stored locally for faster subsequent alignments, and there is an option to annotate the database with ORF positions and frameshift events.

## Installation
### Dependencies
    - biocmanager: 
        - install.packages("BiocManager")
    - ggtree: 
        - BiocManager::install("ggtree")
    - mafft: 
        - sudo apt install mafft
    - raxml: 
        - sudo apt install raxml
    - Edirect:   
        - sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

## Use

    `bash bin/dreamcatch.sh PHYLOGENY_TYPE RUN_TYPE LIST_OF_ACCESSION_NUMBERS (ORF_NAME) (translate)`

    PHYLOGENY_TYPE: 'fast' or 'full'
        fast: quick phylogenetic determination with no bootstrapping
        full: Full phylogenetic determination using 1,000 bootstraps

    RUN_TYPE: 'nucl' or 'amino'
        nucl: Searching NCBI nucleotides database for sequences. Addition of the (translate) option will translate using the standard codon table.
        amino: searches NCBI protein database for sequences.
    
    LIST_OF_ACCESSION_NUMBERS: path to text file containing accession numbers. Included ONLY the accession number per line. Be sure the number matches the intended entry.

    ORF_NAME: Optional argument that will subset the sequence by a START and STOP position annotated in the *_list.csv file. For plain nucleotide sequences and amino acid sequences, the annotation format is `ORF_NAME[START:STOP], ...`. For nucleotide sequences with a frameshift, the format is `ORF_NAME[START:FRAMESHIFTNUMBER_FRAMESHIFPOSITION:STOP], ...". Be sure that the FRAMESHIFTNUMBER is appropriately signed (1 vs -1) and the FRAMESHIFTPOSITION is the position in the original nucleotide sequence, not the subset sequence.

    translate: Optional arugment to convert nucleotide sequence to translate a nucleotide sequence to the amino acid sequence.

### Subsetting sequences

To subset a sequence *without* annotating a specific ORF, add the start and stop positions in the LIST_OF_ACCESSION_NUMBERS file (Example: `EU151723,20,568` will take only positions 20 through 568 for the alignment process).