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
    sequences = []
    with open(fastafile, "r") as f:
        ls = f.readlines()
        for i in ls:
             sequences.append(i.rstrip("\n"))

    seq_dic = {}
    dupl=1
    for i in range(len(sequences)):
        if len(sequences[i]) > 0 and sequences[i][0] == ">" and len(sequences[i+1]) > 0:
            if sequences[i] not in seq_dic.keys():
                seq_id = sequences[i]
            else:
                seq_id = sequences[i]+"_duplicate_"+str(dupl)
                dupl+=1
            seq_seq = sequences[i+1]
            seq_dic[seq_id] = seq_seq
        elif len(sequences[i]) > 0 and sequences[i][0] != ">" and sequences[i-1][0] != ">" :
            seq_seq += sequences[i]
            seq_dic[seq_id] = seq_seq
    if dupl > 1:
        print("Detected multiple instances of sequences with the same name. Only the first instance will be keep the original name.")
        
    return seq_dic

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
