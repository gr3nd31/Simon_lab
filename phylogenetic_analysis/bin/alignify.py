import argparse
seqs_fil = "sequences_aligned.fasta"

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequences", help = "Path to the sequences file")

parser.parse_args()

args = parser.parse_args()
if args.sequences:
    seqs_fil = args.sequences


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

seqs = read_fasta(seqs_fil)
magicNumber = len(seqs[list(seqs.keys())[1]])
pater = "Name,Position,Base\n"

for i in range(0,magicNumber):
    for j in list(seqs.keys()):
        pater+=j[1:]+","+str(i)+","+seqs[j][i]+"\n"

#print(pater)
with open(seqs_fil.replace(".fasta", "_matrix.csv"), 'w') as f:
        f.write(pater)
