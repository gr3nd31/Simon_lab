import argparse

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

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence and structure file")
parser.add_argument("-S", "--StartPosition", help = "Position number to start on.", default=1)
parser.add_argument("-n", "--numberOfAminoAcids", help = "Number of N-terminal amino acids to 'optimize'", default=6)
parser.add_argument("-t", "--typeOfNucleicAcid", help="Which nucleic acid should be optimized (AUGC).", default="A")
parser.add_argument("-m", "--minimizedNucleicAcid", help="Which nucleic acid should be optimized (AUGC).", default="None")
parser.add_argument("-o", "--out", help="Name of output file", default="optiOut.fasta") 
parser.parse_args()
args = parser.parse_args()

runIt = True

# Names the output file
if args.out:
    outFile = args.out

if args.numberOfAminoAcids:
    try:
        nAA=int(args.numberOfAminoAcids)
    except:
        print("Unable to parse given number of amino acids. Please supply an integer.")
        runIt=False

if runIt and args.StartPosition:
    try:
        theStart=int(args.StartPosition)
    except:
        print("Unable to read the start position. Give an integer smaller than the length of the sequence.")


if runIt and args.sequence:
    try:
        theReads = read_fasta(args.sequence)
    except:
        print("Unable to open the target fasta file. Try again?")
        runIt=False

if runIt and args.typeOfNucleicAcid.upper() in ["A", "C", "U", "G", "T"]:
    theMaxer=args.typeOfNucleicAcid.upper().replace("T", "U")
    if args.minimizedNucleicAcid:
        if args.minimizedNucleicAcid.upper() != theMaxer and args.minimizedNucleicAcid.upper() in ["A", "C", "U", "G"]:
            theMiner=args.minimizedNucleicAcid.upper()
        elif args.minimizedNucleicAcid != "None":
            print("Given minimizing nucleic acid not in list of acceptable bases (AUGC). Try again?")
            runIt=False
        else:
            theMiner="None"
else:
    print("Given nucleic acid not in list of acceptable bases (AUGC). Try again?")
    runIt=False

if runIt:
    print("Resetting the first "+str(nAA)+" amino acids to maximize "+theMaxer+"'s and minimize "+theMiner+"'s.")
    theEnd=theStart+(3*nAA)
    for i in theReads.keys():
        theSeq = theReads[i]
        theSeq = theSeq.replace("T", "U")
        print("Optimizing sequence: "+i[1:])
        newString=""
        if theStart>1:
            newString+=theSeq[:theStart]
        for j in range(3,len(theSeq[theStart-1:theEnd-1])+1,3):
            theAmino=codon_table[theSeq[j-3:j]]
            possibleCodons = [j for j in codon_table if codon_table[j] == theAmino]
            if len(possibleCodons)> 1:
                possibleCodons=sorted(possibleCodons, key = lambda s: s.count(theMaxer), reverse=True)
                if args.minimizedNucleicAcid != "None":
                    possibleCodons=sorted(possibleCodons, key = lambda s: s.count(theMiner), reverse=False)
            newString+=possibleCodons[0]          
        newString+=theSeq[(theEnd+theStart)-2:]
        newString=newString.replace("U", "T")
        with open(outFile, 'a') as f:
                f.write(i+"_max"+theMaxer+"_min"+theMiner+"\n"+newString+"\n")