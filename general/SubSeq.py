import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequences", help = "Path to the sequence file.") #
parser.add_argument("-q", "--query", help = "String of the (regex) sequence you want to find") #
parser.add_argument("-u", "--u_output", help = "Output data as 'U'", action='store_true') #
parser.add_argument("-o", "--outputFile", help = "Name of the output file.") #
parser.parse_args()
args = parser.parse_args()

runIt=True
if args.outputFile:
    the_out=args.outputFile
else:
    the_out="SubSeqs.fasta"

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
    theReads=read_fasta(args.sequences)
except:
    print("Unable to open fasta file. Try again?")
    runIt=False

try:
    k=args.query
    k=k.replace("U", "T")
    p = re.compile(k)
except:
    print("Unable to parse given subsequence. Try again?")
    runIt=False

if runIt:
    hits=0
    print("Searching "+str(len(theReads.keys()))+" sequences for: "+k)
    for i in theReads.keys():
        if i == ">CY1":
            print(theReads[i])
        outfile=""
        trick=theReads[i]
        trick=trick.replace("U", "T")
        for m in p.finditer(trick):
            start = m.start()
            end = m.end()
            hits+=1
            outfile+=i+"_"+str(start)+":"+str(end)+"\n"
            tip=theReads[i][start:end]
            if args.u_output:
                tip=tip.replace("T", "U")
            outfile+=tip+"\n"
            with open(the_out, 'a') as f:
                f.write(outfile)
    if hits > 0:
        print("Found sequence "+str(hits)+" times.")
    else:
        print("Sequence not found. Try again?")