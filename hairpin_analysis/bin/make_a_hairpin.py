import argparse
import random
import RNA
import os


def get_pairs(dotBra, rna):
    braDot = dotBra[::-1]

    hitsLib = {'GC':0,
               'AU':0,
               'GU':0,
                'other':[]}
    if len(dotBra) != len(rna):
        print("Sequences are not the same length. Aborting.")
    elif dotBra.count('(') != dotBra.count(')'):
        print("Structure pair is not closed. Aborting.")
    else:
        hits = 0
        for i in range(0,len(dotBra)):
            if dotBra[i] == "(":
                hits += 1
                stih = 0
                for j in range(0,len(braDot)):
                    if braDot[j] == ')':
                        stih += 1

                    if stih == hits:
                        pair = rna[i]+rna[::-1][j]
                        break

                if pair in hitsLib.keys():
                    hitsLib[pair]+=1
                elif pair[::-1] in hitsLib.keys():
                    hitsLib[pair[::-1]]+=1
                else:
                    hitsLib['other'].append(pair)
    return hitsLib

def bulge_count(dotBra):
    b_c = {"apicals": 0,
           "l_bulges": 0,
           "r_bulges": 0}
    b_c["apicals"] = len(re.findall(dotBra, "\\(\\.+\\)"))
    b_c["l_bulges"] = len(re.findall(dotBra, "\\(\\.+\\("))
    b_c["r_bulges"] = len(re.findall(dotBra, "\\)\\.+\\)"))
    return b_c



parser = argparse.ArgumentParser()
parser.add_argument("-n", "--numberOfIterations", help="Number of time to generate the hairpin")
parser.add_argument("-s", "--sequence", help = "Path to the sequence file. If blank, a random hairpin is generated") #
parser.add_argument("-S", "--StartSequence", help="Relative start position of sequence (0-1). Default = 0.") #
parser.add_argument("-b", "--bulgeType", help="Type of bulge (symmetrical, left, right). If blank, no bulges are intentionally created")
parser.add_argument("-B", "--BulgePosition", help="Relative position on the stem for the bulge (0-1). Default 0.5") #
parser.add_argument("-c", "--BulgeSize", help="Number of nucleotides to make the bulge (default 5)") #
parser.add_argument("-p", "--PairedPercent", help = "Percent of hairpin that is paired (0-1, default 1)") #
parser.add_argument("-l", "--length", help="Length of the hairpin (default = 30)") #
parser.add_argument("-a", "--apicalLoop", help = "Length of sequence in the apical loop (default 5)") #
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()
# Assigns the length of the hairpin
if args.length:
    hp_length = int(args.length)
else:
    hp_length = 30

# Pulls a target sequence
if args.sequence:
    seq_file = open(args.sequence, "r")
    sequence = seq_file.read().strip()
    seq_file.close()
    fasta = "hairpin_"+args.sequence
else:
    sequence = ""
    fasta = "hairpin_"+str(hp_length)

# Determines the size of the apical loop
if args.apicalLoop:
    ap_length = int(args.apicalLoop)
else:
    ap_length = 5

# Names the output file
if args.out:
    outFile = args.out
else:
    outFile = "data.csv"

# Determines the percent of stem nt's that are perfectly paired
if args.PairedPercent:
     paired_percent = int(float(args.PairedPercent)*100)
else:
     paired_percent = 100

# Determines the relative position of the target sequence
if args.StartSequence:
    starting=int(float(args.StartSequence)*hp_length)
else:
    starting=int(0.5*(hp_length-len(sequence)))

# Determines the relative position of the bulge
if args.BulgePosition:
    left_bulging=int(float(args.BulgePosition)*hp_length)
    where_the_b=str(args.BulgePosition)
    right_bulging=hp_length-left_bulging
else:
    left_bulging=int(0.5*(hp_length-len(sequence)))
    right_bulging=hp_length-left_bulging
    where_the_b="0.5"

# Determines the size of the bulge
#if symmetrical, the size is the same on both side
#if asymmetrical, the size is only for the one side
if args.BulgeSize:
    bulgeSize = int(args.BulgeSize)
else:
    bulgeSize = 5

nt_set = ["A", "U", "G", "C"]
ops_set = {"A":"U", "U":"A", "C":"G", "G":"C"}

if args.numberOfIterations:
    nnum = int(args.numberOfIterations)+1
else:
    nnum = 1

for iter in range(0, nnum):
# Defines the sequence of the bulge
    leftBulge = ""
    rightBulge = ""
    baseBulge = ""
    twin=-2
    syms=False
    if args.bulgeType:
        codpiece = args.bulgeType
        twin=-2
        while len(baseBulge) < bulgeSize:
            baseBulge+=random.choice(nt_set)

        if args.bulgeType == "left":
            leftBulge = baseBulge
        elif args.bulgeType == "right":
            twin=-2
            rightBulge = baseBulge
        elif args.bulgeType == "symmetrical":
            syms=True
            leftBulge = baseBulge
            for i in leftBulge[::-1]:
                possNTC = random.choice(nt_set)
                while possNTC == ops_set[i]:
                    possNTC = random.choice(nt_set)
                rightBulge += i
        else:
            print("Incorrect bulge type given.")
    else:
        codpiece="None"
        
#--------

    hp_seq = ""
    hidden_seq=""
    seq_added = False
    bulge_added = False
    bulge_ended = False

    counter=0
    counting=False

    while len(hp_seq)-len(leftBulge) < hp_length:
        if len(hp_seq) >= starting and not seq_added:
            seq_added = True
            hp_seq+=sequence
            hidden_seq+=sequence
            if counting:
                counter+=len(sequence)
        elif len(hp_seq) > left_bulging and not bulge_added and len(leftBulge) > 0:
            bulge_added = True
            hp_seq+=leftBulge
            counting=True
        else:
            possNT=random.choice(nt_set)
            hp_seq+=possNT
            hidden_seq+=possNT
            if counting:
                counter+=1

    revc = ""
    timmy=0
    for i in range(len(hidden_seq)-1,twin,-1):
        check = random.randint(0,100)
    # Accounts for differences in left bulges
        if hidden_seq[i] != hp_seq[i+timmy] and not bulge_ended and not syms:
            bulge_ended = True
            revc+=rightBulge
            timmy=1
    
    # Adds the right bulge
        elif len(revc) >= right_bulging and not bulge_ended and not syms:
            bulge_ended = True
            revc+=rightBulge
            timmy=1
    # Adds the right bulge to a symmetrical hairpin
        elif counter== 0 and syms:
            bulge_ended = True
            revc+=rightBulge
            timmy=1
    # Adds the complementary base
        elif check <= paired_percent:
            revc+=ops_set[hidden_seq[i+timmy]]
        else:
            revc += random.choice(nt_set)
        counter-=1

    ap_seq=""
    for i in range(0, ap_length):
        ap_seq+=random.choice(nt_set)

    hp_seq+=ap_seq
    hp_seq+=revc
    if os.path.isfile(outFile):
        data = ""
    else:
        data = "Name,Complementarity,ApicalSize,ApicalSeq,Bulge,BulgeSize,BulgePosition,BulgeSeq,StemLength,Sequence,Structure,Length,bp,GC,dG,dG_Length,PE,APE,GC_pairs,AU_pairs,GU_pairs,GC_pair_percent,AU_pair_percent,GU_pair_percent,apicals,left_bulges,right_bulges\n"

    fc = RNA.fold_compound(hp_seq)
    fc.pf()
    
    pairs = get_pairs(fc.mfe()[0], hp_seq)
    bulge_counts = bulge_count(fc.mfe()[0])


    data+=fasta+"," #Adds a general name
    data+=str(paired_percent)+"," #Adds complementarity score
    data+=str(ap_length)+"," # Adds the size of the apical loop
    data+=ap_seq+"," # Adds the apical sequence
    data+=codpiece+"," #Adds the bulge type
    data+=str(len(baseBulge))+"," #Adds Bulge size
    data+=where_the_b+"," #Add Bugle position
    data+=baseBulge+"," #Adds bulge sequence
    data+=str(hp_length)+"," #Adds the length of hairpin stem
    data+=hp_seq+"," #Adds the primary sequence
    data+=fc.mfe()[0]+"," # Adds the structure sequence in dot-bracket
    data+=str(len(hp_seq))+"," # Adds the total length of the sequence
    data+=str(fc.mfe()[0].count("("))+"," # Adds the number of basepaired bases
    data+=str(round((hp_seq.count("G")+hp_seq.count("C"))/len(hp_seq),2))+"," # Adds teh GC content
    data+=str(round(fc.mfe()[1],3))+"," #Adds the MFE
    data+=str(round(fc.mfe()[1]/len(hp_seq),2))+"," # Adds the MFE/length ratio
    data+=str(fc.positional_entropy()).replace(",","")+"," #Adds all the PEs
    trick=list(fc.positional_entropy())
    data+=str(sum(trick[1:])/(len(trick)-1))+"," #Adds APE
    data+=str(pairs['GC'])+"," #Adds number of GC pairings
    data+=str(pairs['AU'])+"," #Adds number of AU pairings
    data+=str(pairs['GU'])+"," #Adds number of GU pairings
    if fc.mfe()[0].count("(") > 0:
        pair_count = fc.mfe()[0].count("(")
    else:
        pair_count = 1
    data+=str(pairs['GC']/pair_count)+"," #Adds percent of GC pairings
    data+=str(pairs['AU']/pair_count)+"," #Adds percent of AU pairings
    data+=str(pairs['GU']/pair_count)+"," #Adds percent of GU pairings
    data+=str(bulge_counts["apicals"])+"," #Adds the number of apical loops
    data+=str(bulge_counts["l_bulges"])+"," #Adds the number of 5 prime bulges
    data+=str(bulge_counts["r_bulges"])+"\n" #Adds the number of 3 prime bulges

    with open(outFile, 'a') as f:
        f.write(data)