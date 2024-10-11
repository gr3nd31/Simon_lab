import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file.") #
parser.add_argument('-p', '--percentDifferent', help="Percent (0-1) of possible differences that should be made. Default = 0.5")
parser.add_argument('-c' '--coding', help= "Keep coding sequence the same. Should be an integer indicating the position where coding begins. Default is 0 (first base)")
parser.add_argument('-e', '--endCoding', help= "Position to end the coding on.")
parser.add_argument('-u', '--usageTable', help = 'Path to codon usage table')
parser.add_argument('-n', '--numberOfIterations', help="Number of iterations to make. Default is 100.")
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()
runIt = True
caring = False

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

nt_table = ["U", "G", "C", "A"]

# Pulls a target sequence
if args.sequence:
    try:
        seq_file = open(args.sequence, "r")
        sequence = seq_file.read().strip()
        sequence = sequence.upper().replace(" ", "")
        sequence = sequence.replace("T", "U")
        seq_file.close()
        fasta = args.sequence
    except:
        print("Unable to read sequence file. Abortin.")
        sequence = ""
        fasta = str(hp_length)
        runIt = False
else:
    sequence = ""
    fasta = "hairpin_"+str(hp_length)
    runIt = False

if args.usageTable:
    try:
        seq_file = open(args.usageTable, "r")
        cTable = seq_file.read().strip()
        cTable = cTable.replace("T", "U"),
        tabs = cTable[0].splitlines()
        usageTable = {}
        for i in tabs:
            usageTable[i.split(" ")[0]] = float(i.split(" ")[1])
        caring = True
        seq_file.close()
    except:
        print("Unable to parse codon table. Mutating without caring.")

# Names the output file
if args.out:
    outFile = args.out
else:
    outFile = "Mutated_"+args.sequence

if args.c__coding:
    keeping = True
    try:
        frame_number = int(args.c__coding)-1
    except:
        print("Integer position not parseable. Beginning coding from the first base.")
        frame_number = 1
else:
    keeping = False
    frame_number = 0

if args.endCoding and keeping:
    try:
        final_number = int(args.endCoding)
    except:
        print("Final position not parseable. Ending coding at the last codon.")
        final_number = len(sequence)
else:
    final_number = len(sequence)

if args.numberOfIterations:
    try:
        numberOfIterations = int(args.numberOfIterations)
    except:
        print("Unable to parse iterations number, defaulting to 100")
        numberOfIterations = 100
else:
    numberOfIterations = 100

if args.percentDifferent:
    try:
        percentDiff = float(args.percentDifferent)
    except:
        print("Unable to parse given percent difference, defaulting to 50% (0.5)")
        percentDiff = 0.5
else:
    percentDiff = 0.5

up_seq=sequence[0:frame_number]
down_seq=sequence[final_number:len(sequence)]

if keeping:
    new_seq = ""
    for i in range(0,len(sequence), 3):
        target_seq = sequence[frame_number:final_number][i:i+3]
        if len(target_seq) == 3:
            new_seq+=codon_table[target_seq]
        else:
            down_seq=target_seq+down_seq
    the_seq = new_seq
else:
    the_seq = sequence


if runIt:
    finalOut = ">"+fasta.replace(".txt", "")+"_"+"Original\n"+sequence+"\n"
    # Calculate the number of changes to be made
    requestedChanges = round(percentDiff*len(the_seq))
    if requestedChanges == 0:
        print("Requested changes is 0. Moving to 1.")
        requestedChanges+=1
    else:
        print("Making "+str(requestedChanges)+" changes...")

    # Iterate through the number of iterations requested
    for i in range(0, numberOfIterations):
        iter_seq = the_seq

        # Append a sequence name
        finalOut+=">"+fasta.replace(".txt", "")+"_"+str(requestedChanges)+"-changes"+"_"+str(i+1)+"\n"

        # Mutates codons
        if keeping:
            attempts = 0
            changes = []
            while len(changes) < requestedChanges and attempts <= 100:
                hit = random.sample(range(0, len(iter_seq)), 1)[0]
                if hit*3 not in changes:
                    possibleCodons = [j for j in codon_table if codon_table[j] == iter_seq[hit]]
                    if len(possibleCodons) > 1:
                        changes.append(hit*3)
                    else:
                        if attempts == 99:
                            print("Passed 100 attempts to mutate hairpin with coding. Try decreasing the requested number of changes.")
                        attempts+=1
           
            finalOut+=up_seq
            stim = sequence[len(up_seq):len(sequence)-len(down_seq)]
            for k in range(0,len(stim),3):
                if k in changes:
                    possibleCodons = [j for j in codon_table if codon_table[j] == iter_seq[int(k/3)]]
                    # sorts codons by usage
                    if caring:
                        usage = []
                        summation = 0
                        for u in possibleCodons:
                            summation += usageTable[u]
                        for u in possibleCodons:
                            usage.append((u, usageTable[u]/summation))
                        usage.sort(key=lambda tup: tup[1])
                        tack = 0
                        for u in range(0,len(usage)):
                            the_hit = usage[u][0]
                            the_num = usage[u][1]+tack
                            usage[u]=(the_hit, int(the_num*100))
                            tack=the_num
                    trip = stim[k:k+3]
                    while trip == stim[k:k+3]:
                        if caring and len(possibleCodons) > 2:
                            the_bean = random.randint(1,100)
                            for u in usage:
                                if u[1] > the_bean:
                                    trip = u[0]
                                    break
                        else:
                            trip = random.choice(possibleCodons)
                    finalOut+=trip
                else:
                    finalOut+=stim[k:k+3]
            finalOut+=down_seq+"\n"

        # Mutates bases
        else:
            changes = random.sample(range(0, len(sequence)), requestedChanges)
            for i in changes:
                hit = iter_seq[i]
                new = iter_seq[i]
                while new == hit:
                    new = random.choice(nt_table)
                stim = list(iter_seq)
                stim[i] = new
                iter_seq = "".join(stim)
            finalOut+=iter_seq+"\n"
with open(outFile, 'a') as f:
    f.write(finalOut)