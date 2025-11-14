import argparse

seqs_fil = "sequences.fasta"
codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequences", help = "Path to the sequences file")
parser.add_argument("-t", "--tag", help = "Tag to add to sequence name")
parser.add_argument("-o", "--outFile", help = "Name of the output file. By default, it is the in-file.")
parser.parse_args()

args = parser.parse_args()
if args.sequences:
    seqs_fil = args.sequences

if args.tag:
    id_tag = "_"+args.tag
else:
    id_tag = ""

if args.outFile:
    outFile = args.outFile
else:
    outFile = args.sequences


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

seqs = read_fasta(seqs_fil)
new_set = ""
for i in seqs:
    #print(seqs[i])
    new_set=new_set+i+id_tag+"\n"
    nner = seqs[i].upper()
    if "U" in nner.upper():
        nner = nner.replace("U", "T")
    end = len(nner) - (len(nner) %3) - 1
    for j in range(0,end, 3):
        codon = nner[j:j+3]
        if codon in codon_table:
            new_set += codon_table[codon]
        else:
            new_set += "N"
    new_set+="\n"

#print(new_set)
with open(outFile, 'w') as f:
        f.write(new_set)
