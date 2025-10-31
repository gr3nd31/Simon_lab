import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the fasta file.")
parser.add_argument("-n", "--number_to_sample", help = "Number of reads to be sampled.", default=1000000)
parser.add_argument("-o", "--out", help="Name of output file", default="sampled_reads.fasta")
parser.parse_args()
args = parser.parse_args()

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

runIt=True
if args.sequence:
    try:
        fastas=read_fasta(args.sequence)
    except:
        print("Unable to read sequence file. Remember to give a fasta file.")
        runIt=False
else:
    print("No sequence file given. Remember to give a fasta file.")
    runIt=False

if args.number_to_sample:
    try:
        n=int(args.number_to_sample)
    except:
        print("Unable to read given sample number. Remember to give an integer.")

if args.out:
    outFile=args.out

print("Randomly sampling "+str(n)+" reads.")
if len(fastas.keys()) < n:
    print("Requested sampling number is greater than the number of reads in the fasta ("+str(len(fastas.keys()))+"). Try using a smaller integer.")
    runIt=False

if runIt:
    idList=[]
    kList=list(fastas.keys())
    while len(idList) < n:
        x=random.choice(kList)
        if x not in idList:
            idList.append(x)
            kList.remove(x)
    for i in idList:
        outList=i+"\n"+fastas[i]+"\n"
        with open(outFile, 'a') as f:
            f.write(outList)
