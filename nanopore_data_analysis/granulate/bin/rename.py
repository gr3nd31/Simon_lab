import argparse
# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help = "Path to input file")
parser.add_argument("-o", "--out", help = "Out file. Defulat is input file with 'new_' added")
parser.parse_args()
args = parser.parse_args()
runIt=True

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
    x=read_fasta(args.input)
except:
    print("Unable to open reads file. Aborting.")
    runIt=False

if args.out:
    outFile = args.out
else:
    outFile="new_"+args.input

if runIt:
    try:
        for i in x:
            outstring=i.split(" ")[0]+"\n"+x[i]+"\n"
            with open(outFile, "a") as f:
                f.write(outstring)
    except:
        print("Unable to write file using "+outFile+" file name. Aborting")
