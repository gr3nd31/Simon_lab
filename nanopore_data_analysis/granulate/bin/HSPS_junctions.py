import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequences", help = "Path to the fasta file") #
parser.add_argument("-j", "--junctionCSV", help="Path to the junction CSV file with the junctions to be sliced") #
parser.add_argument("-o", "--out", help="Name of output file", default="output.csv") #
parser.parse_args()
args = parser.parse_args()


def read_fasta(fastafile):
    """
    Reads a fasta file and returns a dictionary with sequence
    number as keys and sequence code as values
    """
    sequences = {}
    with open(fastafile, "r") as f:
        ls = f.read()
    ls.rstrip("\n")
    split_reads=ls.split(">")
    for i in split_reads:
        j=i.split("\n")
        if j[0] != "":
            seqName=j[0]
            theSeq=""
            for k in j[1:]:
                theSeq+=k
            if seqName not in sequences.keys():
                sequences[seqName]=theSeq
            else:
                print("Duplicate seqID found for: "+seqName[1:])
    return sequences

try:
    sequences = read_fasta(args.sequences)
except FileNotFoundError:
    print("Error: Fasta file not found")
    exit(1)

try:
    junctions = pd.read_csv(args.junctionCSV)
except FileNotFoundError:
    print("Error: Junction CSV file not found")
    exit(1)

try:
    junctions['Sequence'] = junctions.apply(lambda row: sequences[row['readId']][min([row['gapStart'], row['gapEnd']]):max([row['gapStart'], row['gapEnd']])-1], axis=1)
except KeyError:
    print("Error: Some read IDs in the junction CSV file are not found in the fasta file")
    exit(1)

junctions.to_csv(args.out, index=False)