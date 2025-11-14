import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file.") #
parser.add_argument("-w", "--windowSize", help = "Length of the window size.") #
parser.add_argument("-o", "--outputFile", help = "Name of the output file.") #
parser.parse_args()
args = parser.parse_args()

if args.windowSize:
    window_size = int(args.windowSize)
else:
    window_size = 21

if args.outputFile:
    the_out=args.outputFile
else:
    the_out="charge_data.csv"

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

charge_tabel = {"R":1, "H":1, "K":1,
                "D":-1, "E":-1,
                "S":0, "T":0, "N":0, "Q":0,
                "C":0, "G":0, "P":0,
                "A":0, "V":0, "I":0, "L":0,
                "M":0, "F":0, "Y":0, "W":0}

# Pulls a target sequence
if args.sequence:
    try:
        seq_file = read_fasta(args.sequence)
    except:
        print("Unable to read sequence file. Abortin.")
        runIt = False

the_csv="Chain,Position,Charge,Residues\n"
for i in seq_file:
    #print(seq_file[i])
    the_seq = seq_file[i].upper()
    for j in range(0,len(seq_file[i])):
        window_charge = 0
        icky = 0
        for x in range(j,j+window_size):
            try:
                #print(the_seq[x])
                window_charge+=charge_tabel[the_seq[x]]
                icky+=1
            except:
                pass
        #print("window "+str(j)+" average charge is "+str(round(window_charge/window_size, 4)))
        the_csv+=i.replace(">", "")+","+str(j)+","+str(round(window_charge/window_size, 4))+","+str(icky)+"\n"
        #print(icky
#print(the_csv)
with open(the_out, 'w') as f:
    f.write(the_csv)