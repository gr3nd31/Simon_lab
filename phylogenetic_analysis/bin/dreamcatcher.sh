#! /bin/bash

#$1 = Type of run: `fast` vs `full`
#echo $1

#$2 = Type of analysis: `nucl` for nucleotide or `amino` for amino acids (NOT FUNCTIONAL YET)
#echo $2

#$3 = Path to list file
#echo $3

#$4 = Name of ORF to compare (MUST BE ANNOTATED IN THE DB FILE AS `ORFNAME[START:STOP],...`)
#echo $4

#$5 = translate


if [ "$2" == "nucl" ]; then
    run_it=true
    echo "Running nucleic acid alignment"
    # Creates sequence repository and id list
    if [ ! -d "nucleotide_seqs" ]; then
        echo "Nucleotide repository not detected. Creating directory"
        mkdir nucleotide_seqs
        echo "Accession_id,Name,ORFs,Tag" > nucleotide_list.csv
    fi

    # Checks that nucleotide list exists and makes it if it doesn't
    if [ ! -f "nucleotide_list.csv" ]; then
        echo "Accession_id,Name,ORFs,Tag" > nucleotide_list.csv
    fi

    # Looks for missing accession numbers
    # Yes, this can be done in bash, but I'm not good enough to split by commas and since I'm using R anyway...
    Rscript bin/find_missing.R $3 nucleotide_list.csv

    # Checks if sequence is already present
    if [ -f pull_list.txt ]; then
        # Splits pull list and returns a unqiue list
        array=($(echo $(cat pull_list.txt) | tr "," "\n"))
        uniqs_arr=($(for acc in "${array[@]}"; do echo "${acc}"; done | sort -u))
        echo "Beginning sequence fetching."
        for word in "${uniqs_arr[@]}"; do
        # If not present, sequence is downloaded and added to sequence logs
        if [ ! -f nucleotide_seqs/$word ]; then
            echo "$word sequence not found. Retrieving..."
            esearch -db nucleotide -query $word | efetch -format fasta > nucleotide_seqs/$word
            bim=$word
            bim+=",,NONE,NONE"
            echo $bim >> nucleotide_list.csv
        fi
    done
    echo "Sequence fetch complete."
    rm pull_list.txt
    fi

    # Generates alignment fasta
    Rscript bin/find_names_and_frames.R $3 nucleotide_list.csv $4 $5

    if [ "$5" == "translate" ]; then
        echo "Converting nucleotides to amino acids..."
        python3 bin/nucl_to_aa.py -s sequences.fasta
        translated=true
    fi

    # Runs mafft alignment and removes the raw sequence file
    echo "Beginning alignment"
    mafft --auto sequences.fasta > sequences_aligned.fasta
    #rm sequences.fasta

    # Generates a phylogenetic tree
    echo "Generate relationships"
    if [ -f "RAxML_bestTree.sequences_tree" ]; then
        rm RAxML*
    fi

    if [ "$5" == "translate" ]; then
        if [ "$1" == "full" ]; then
            raxmlHPC -T 8 -m PROTGAMMAAUTO -p 144 -f a -x 144 -# 1000 -s sequences_aligned.fasta -n sequences_tree 
        else
            raxmlHPC -s sequences_aligned.fasta -n sequences_tree -m PROTGAMMAAUTO -p 12
        fi
    else
        if [ "$1" == "full" ]; then
            raxmlHPC -T 8 -m GTRGAMMAI -p 144 -f a -x 144 -# 1000 -s sequences_aligned.fasta -n sequences_tree 
        else
            raxmlHPC -d -p 12 -m GTRGAMMAI -s sequences_aligned.fasta -n sequences_tree
        fi
    fi
#------------------------------------------------------------------
elif [ "$2" == "amino" ]; then
    run_it=true
    echo "Running amino acid alignment"
    # Creates sequence repository and id list
    if [ ! -d "aa_seqs" ]; then
        echo "Amino acid repository not detected. Creating directory"
        mkdir aa_seqs
        echo "Accession_id,Name,ORFs,Tag" > amino_acid_list.csv
    fi

    # Creates sequence log if missing
    if [ ! -f "./amino_acid_list.csv" ]; then
        echo "Accession_id,Name" > amino_acid_list.csv
        for fasta_file in $(ls aa_seqs/); do
            bim=$(echo $fasta_file | cut -d "/" -f 2)
            bim+=","
            echo $bim >> amino_acid_list.csv
        done
    fi
    # Looks for missing accession numbers
    # Yes, this can be done in bash, but I'm not good enough to split by commas and since I'm using R anyway...
    Rscript bin/find_missing.R $3 amino_acid_list.csv

    # Checks if sequence is already present
    if [ -f pull_list.txt ]; then
        # Splits pull list and returns a unqiue list
        array=($(echo $(cat pull_list.txt) | tr "," "\n"))
        uniqs_arr=($(for acc in "${array[@]}"; do echo "${acc}"; done | sort -u))
        echo "Beginning sequence fetching."
        for word in "${uniqs_arr[@]}"; do
        # If not present, sequence is downloaded and added to sequence logs
        if [ ! -f aa_seqs/$word ]; then
            echo "$word sequence not found. Retrieving..."
            esearch -db protein -query $word | efetch -format fasta > aa_seqs/$word
            bim=$word
            bim+=",,NONE,NONE"
            echo $bim >> amino_acid_list.csv
        fi
    done
    echo "Sequence fetch complete."
    rm pull_list.txt
    fi

    # Generates alignment fasta
    Rscript bin/find_names_and_frames.R $3 amino_acid_list.csv $4

    # Runs mafft alignment and removes the raw sequence file
    echo "Beginning alignment"
    mafft --auto sequences.fasta > sequences_aligned.fasta
    #rm sequences.fasta

    # Generates a phylogenetic tree
    echo "Generate relationships"
    if [ -f "RAxML_bestTree.sequences_tree" ]; then
        rm RAxML*
    fi

    if [ "$1" == "full" ]; then
        raxmlHPC -T 8 -m PROTGAMMAAUTO -p 144 -f a -x 144 -# 1000 -s sequences_aligned.fasta -n sequences_tree 
    else
        raxmlHPC -s sequences_aligned.fasta -n sequences_tree -m PROTGAMMAAUTO -p 12
    fi

#----------------------------------------------------------------
else
    echo "Incorrect alignment type given. Try again please."
    run_it=false
fi

if $run_it; then
    # Creates a tree pdf
    echo "Creating tree"
    Rscript bin/ggtree.R RAxML_bestTree.sequences_tree
    echo Completed
fi