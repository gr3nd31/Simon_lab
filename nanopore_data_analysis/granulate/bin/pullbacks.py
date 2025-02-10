import argparse

default_file = "list.txt"
default_fasta_reads = "all_reads.fasta"
default_out = "out_reads.fasta"
exclude = False

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-l", "--list_to_pull", help = "Text file of read IDs to be pulled")
parser.add_argument("-r", "--reads", help = "Fasta file of all reads")
parser.add_argument("-o", "--out", help = "Text file of reads to be pulled")
parser.add_argument("-e", "--exclude", help = "Exclude these reads instead of pulling them")
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


fastq_reads = read_fasta(default_fasta_reads)

hits = 0
for i in x:
    hit_reads = ">"+i+"\n"
    if ">"+i in fastq_reads.keys():
        hits+=1
        hit_reads+=fastq_reads[">"+i]+"\n"
        with open(default_out, "a") as f:
            f.write(hit_reads)

print("\nNumber of reads requested: " + str(len(x)))
print("Number of reads searched: "+str(len(fastq_reads)))
trick = hits/len(fastq_reads)
print("Number of reads found: "+str(hits)+" ("+str(round(100*hits/len(x)))+"%)")
