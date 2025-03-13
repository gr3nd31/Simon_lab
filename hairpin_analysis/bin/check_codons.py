import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file.") #
parser.add_argument('-c', '--coding', help= "Keep coding sequence the same. Should be an integer indicating the position where coding begins. Default is 0 (first base)")
parser.add_argument('-e', '--endCoding', help= "Position to end the coding on.")
parser.add_argument('-u', '--usageTable', help = 'Path to codon usage table')
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()
runIt = True

codon_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
        }

amino_table = {"I":0, "M":0, "T":0, "N":0, "K":0, "S":0,
               "R":0, "L":0, "P":0, "H":0, "Q":0, "V":0,
               "A":0, "D":0, "E":0, "G":0, "F":0, "Y":0,
               "C":0, "_":0, "W":0}

if args.sequence:
    try:
        seq_file = open(args.sequence, "r")
        sequence = seq_file.read().strip()
        seq_file.close()
        fasta = sequence
    except:
        runIt = False

if args.usageTable:
    try:
        seq_file = open(args.usageTable, "r")
        sequence = seq_file.read().strip()
        seq_file.close()
        use_table = sequence.split("\n")
        catch = {}
        for i in use_table:
            codon = i.split(" ")[0]
            percent = float(i.split(" ")[1])
            catch[codon] = [codon_table[codon], percent, 0, 0]
        
    except:
        runIt = False

# Names the output file
if args.out:
    outFile = args.out
else:
    outFile = "codonUsage_"+args.sequence.split(".")[0]+".csv"

if args.coding:
    try:
        startPoint = int(args.coding)
    except:
        print("Unable to parse given start position. Defaulting to first position")
        startPoint = 0
else:
    startPoint = 0

if args.endCoding:
    try:
        endPoint = int(args.endCoding)
    except:
        print("Unable to parse given start position. Defaulting to first position")
        endPoint = len(fasta)
else:
    endPoint = len(fasta)

if runIt:
    fasta_name = fasta.split()[0]
    fasta_seq = fasta.split()[1]
    fasta = fasta_seq[startPoint:endPoint]
for i in range(0, len(fasta), 3):
        tick = fasta[i:i+3]
        catch[tick][2]+=1
        amino_table[catch[tick][0]]+=1
    
for i in catch:
    if amino_table[catch[i][0]] > 0:
        catch[i][3] = round(catch[i][2]/amino_table[catch[i][0]], 2)
    else:
        catch[i][3] = 0

outString = "codon,amino,count,aCount,expPercent,realPercent\n"
for i in catch:
    outString+=i+","+catch[i][0]+","+str(catch[i][2])+","+str(amino_table[catch[i][0]])+","+str(catch[i][1])+","+str(catch[i][3])+"\n"

with open(outFile, 'w') as f:
    f.write(outString)