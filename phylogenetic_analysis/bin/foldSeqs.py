import argparse
import RNA
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file. If blank, a random hairpin is generated") #
parser.add_argument("-o", "--out", help="Name of output file", default="sequences_folded.fasta") #
parser.parse_args()
args = parser.parse_args()
run_it=True

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
            seqName=">"+j[0]
            theSeq=""
            for k in j[1:]:
                theSeq+=k
            if seqName not in sequences.keys():
                sequences[seqName]=theSeq
            else:
                print("Duplicate seqID found for: "+seqName[1:])
    return sequences

try:
    sequences = read_fasta(args.sequence)
except:
    print("Unable to read fasta file. Please check the path and try again.")
    run_it=False

if run_it:
    for seqID in sequences.keys():
        seq = sequences[seqID].replace("T", "U").upper()
        if re.search(r"[^ACGU]", seq):
            print("Sequence "+seqID[1:]+" contains invalid characters. Please check the sequence and try again.")
        else:
            fc = RNA.fold_compound(seq)
            fc.pf()
            structure = fc.mfe()[0].replace(",", ".")
            #structure = re.sub("\\(\\.+\\)", "A", structure)
            structure = structure.replace("(", "C")
            structure = structure.replace(")", "G")
            structure = structure.replace(".", "U")
            with open(args.out, 'a') as f:
                f.write(seqID+"\n")
                f.write(structure+"\n")
            