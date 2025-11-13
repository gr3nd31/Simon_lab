import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file.") #
parser.add_argument("-S", "--SubSequence", help = "String of the (regex) sequence you want to find") #
parser.add_argument("-e", "--extraSequenceLength", help = "How much extra sequence is okay (assuming you're not great with regex)", default=20)
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
    sequences = []
    with open(fastafile, "r") as f:
        ls = f.readlines()
        for i in ls:
             sequences.append(i.rstrip("\n"))

    seq_id = []
    for i in sequences:
        if len(i) > 0:
            if i[0] == ">":
                seq_id.append(i)

    seq_id_index = []
    for i in range(len(seq_id)):
        seq_id_index.append(sequences.index(seq_id[i]))

    seq_dic = {}
    for i in range(len(seq_id_index)):
        if i == (len(seq_id_index) - 1):
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:]
        else:
            seq_dic[seq_id[i]] = sequences[seq_id_index[i]+1:seq_id_index[i+1]]

    seq_dic_2 = {}
    for keys, values in seq_dic.items():
        seq_dic_2[keys] = "".join(values)

    return seq_dic_2

try:
    theReads=read_fasta(args.sequence)
except:
    print("Unable to open fasta file. Try again?")
    runIt=False

try:
    k=args.SubSequence
    k=k.replace("U", "T")
    p = re.compile(k)
except:
    print("Unable to parse given subsequence. Try again?")
    runIt=False

if runIt:
    hits=0
    print("Searching "+str(len(theReads.keys()))+" sequences for: "+k)
    for i in theReads.keys():
        outfile=""
        trick=theReads[i]
        trick=trick.replace("U", "T")
        for m in p.finditer(trick):
            start = m.start()
            end = m.end()
            if len(k)+args.extraSequenceLength >= end-start:
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