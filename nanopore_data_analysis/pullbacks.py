import argparse

default_file = "list.txt"
default_fasta_reads = "all_reads.fasta"
default_out = "cy1_di"
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

gin = open(default_fasta_reads)
fastq_reads = gin.read()
gin.close()
fastq_reads = fastq_reads.split("\n>")
hits = 0
hit_reads = ""
#print(fastq_reads[0])
#print(x[0])
for j in fastq_reads:
    label = j.split("\n")[0]
    if len(label.split(" ")) > 1:
        label = label.split(" ")[0]
    if ">" in label:
        label = label.split(">")[1]
        #print(label)
    if label in x and not exclude:
        hits +=1
        #print("Found: "+label)
        try:
            hit_reads+=">"+label+"\n"+j.split("\n")[1]+"\n"
        except:
            continue
            print(j)
    elif label not in x and exclude:
        hits +=1
        #print("Found: "+label)
        try:
            hit_reads+=">"+label+"\n"+j.split("\n")[1]+"\n"
        except:
            continue
            print(j)
print("\nNumber of reads requested: " + str(len(x)))
print("Number of reads searched: "+str(len(fastq_reads)))
trick = hits/len(fastq_reads)
print("Number of reads found: "+str(hits)+" ("+str(round(100*hits/len(x)))+"%)")

#print(hit_reads)
with open(default_out, "w") as f:
    f.write(hit_reads)
