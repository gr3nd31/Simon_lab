import argparse
import numpy as np

default_file = "list.txt"
default_fasta_reads = "all_reads.fasta"
default_out = "out_reads.fasta"
exclude = False

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-l", "--list_to_pull", help = "Text file of read IDs to be pulled")
parser.add_argument("-r", "--reads", help = "Fasta file of all reads")
parser.add_argument("-o", "--out", help = "Text file of reads to be pulled")
parser.add_argument("-e", "--exclude", help = "Exclude these reads instead of pulling them", action = "store_true")
parser.parse_args()
args = parser.parse_args()

if args.list_to_pull:
    default_file = args.list_to_pull
    #print("yeah")
if args.reads:
    #print("yeah")
    default_fasta_reads = args.reads
    #print(default_fasta_reads)
if args.exclude:
    exclude=True
if args.out:
    if not args.out.endswith(".fasta"):
        default_out = args.out+".fasta"
    else:
        default_out = args.out

if not exclude:
    print("\nPulling reads listed in: "+default_file)
else:
    print("\nRemoving reads listed in: "+default_file)
print("Pulling reads from: "+default_fasta_reads)

f = open(default_file)
x = f.read()
f.close()
x = x.strip().split("\n")

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


fastq_reads = read_fasta(default_fasta_reads)

hits = 0
if not exclude:
    for i in x:
        hit_reads = ">"+i+"\n"
        if ">"+i in fastq_reads.keys():
            hits+=1
            hit_reads+=fastq_reads[">"+i]+"\n"
            with open(default_out, "a") as f:
                f.write(hit_reads)
elif exclude:
    trick = np.array(x)
    trick=">"+trick
    theDiff = set(fastq_reads.keys())-set(trick)
    for i in theDiff:
        hit_reads = ">"+i+"\n"
        if i in fastq_reads.keys():
            hits+=1
            hit_reads+=fastq_reads[i]+"\n"
            with open(default_out, "a") as f:
                f.write(hit_reads)


if not exclude:
    print("\nNumber of reads requested: " + str(len(x)))
else:
    print("\nNumber of reads requested: " + str(len(fastq_reads)-len(x)))
print("Number of reads searched: "+str(len(fastq_reads)))
trick = hits/len(fastq_reads)
print("Number of reads found: "+str(hits)+" ("+str(round(100*hits/len(x)))+"%)")
