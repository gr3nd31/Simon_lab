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
    sequences = []
    with open(fastafile, "r") as f:
        ls = f.readlines()
        for i in ls:
             sequences.append(i.rstrip("\n"))

    seq_id = []
    for i in sequences:
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

seqs = read_fasta(seqs_fil)
magicNumber = len(seqs[list(seqs.keys())[1]])
pater = "Name,Position,Base\n"

for i in range(0,magicNumber):
    for j in list(seqs.keys()):
        pater+=j[1:]+","+str(i)+","+seqs[j][i]+"\n"

#print(pater)
with open(seqs_fil.replace(".fasta", "_matrix.csv"), 'w') as f:
        f.write(pater)
