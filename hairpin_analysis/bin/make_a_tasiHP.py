import argparse
import random
import RNA
import os
import re

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

def bulge_count(dotBra, rna):
    b_c = {"apicals": 0,
           "apical_local": "",
           "l_bulges": 0,
           "lB_local": "",
           "r_bulges": 0,
           "rB_local": ""}
    b_c["apicals"] = len(re.findall("\\(\\.+\\)", dotBra))
    ap=re.compile("\\(\\.+\\)")
    for found in ap.finditer(dotBra):
        b_c["apical_local"]+=rna[found.start()+1:found.end()-1]+";"
    b_c["l_bulges"] = len(re.findall("\\(\\.+\\(", dotBra))
    ap=re.compile("\\(\\.+\\(")
    for found in ap.finditer(dotBra):
        b_c["lB_local"]+=rna[found.start()+1:found.end()-1]+";"
    b_c["r_bulges"] = len(re.findall("\\)\\.+\\)", dotBra))
    ap=re.compile("\\)\\.+\\)")
    for found in ap.finditer(dotBra):
        b_c["rB_local"]+=rna[found.start()+1:found.end()-1]+";"
    return b_c

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--numberOfIterations", help="Number of time to generate the hairpin", default=1)
parser.add_argument("-s", "--sequence", help = "Path to the sequence file.") #
parser.add_argument("-A", "--Anchor", help='Whether the miRNA targeting site should be anchored to the apical loop (A) or the base of the hairpin (B).', 
                    default="A", choices=['A', 'B'])
parser.add_argument("-a", "--apicalSize", help = "Length of the apical loop.", default=5) #
parser.add_argument("-R", "--Repeats", help="Flag repeats. Give an integer length to check if sequence is repeated.")
parser.add_argument("-D", "--Details", help="Report sequences of apical loops and bulges rather than number.",  action='store_true')
parser.add_argument("-p", "--PairedPercent", help = "Percent of hairpin that is paired (0-1, default 0.6)") #
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()
lcounter = 0
runIt=True

# Option to check is a sequence is repeated in the hairpin
if args.Repeats:
    try:
        repWindow=int(args.Repeats)
        foundRepeat=False
    except:
        print("Unable to parse given integer. Defaulting to 5.")
        repWindow=5
        foundRepeat=False
else:
    repWindow=0
    foundRepeat=False

# Pulls a target sequence
if args.sequence:
    try:
        fastas=read_fasta(args.sequence)
    except:
        print("Unable to read sequence file. Remember to give a fasta file.")
        runIt=False
else:
    print("No sequence file given. Remember to give a fasta file.")
    runIt=False

# Names the output file
if args.out:
    outFile = args.out
else:
    outFile = "data.csv"

# Defines the size of the apical loop
if args.apicalSize:
    try:
        apicalSize = int(args.apicalSize)
    except:
        print("Unable to parse the apical loop size. Enter an integer. Aborting")
        runIt=False

# Determines the percent of stem nt's that are perfectly paired
if args.PairedPercent:
     paired_percent = int(float(args.PairedPercent)*100)
else:
     paired_percent = 60

nt_set = ["A", "U", "G", "C"]
ops_set = {"A":"U", "U":"A", "C":"G", "G":"C"}

if args.numberOfIterations:
    nnum = int(args.numberOfIterations)+1
else:
    nnum = 1

if runIt:
    if args.Anchor == "A":
        print("Anchoring sequence to the apical loop.")
        aType="apical"
    else:
        print("Anchoring sequence to the base of the hairpin.")
        aType="basal"
    for i in fastas.keys():
        print("Generating hairpins of sequence: "+i[1:])
        fullSeq=fastas[i]
        fullSeq=fullSeq.upper()
        if "T" in fullSeq:
            fullSeq=fullSeq.replace("T", "U")

        if (len(fullSeq)-33)%21 != 0:
            pause=input("Given sequence doesn't match the 22+11+21*n format of a tasiRNA sequence. Continue (Y/n)? ")
            if pause.upper() == "N":
                print("Skipping this sequence.")
                continue

        for iter in range(0, nnum-1):
            altSeq=""
            lcounter+=1
            if lcounter%100 == 0:
                print("Folding read: "+str(lcounter))

            if args.Anchor == "A":
                stemSeq=fullSeq[apicalSize:]
                apicalSeq=fullSeq[:apicalSize]
                for j in range(1,len(stemSeq)+1):
                    check = random.randint(0,100)
                    if check <= paired_percent or j < 4:
                        altSeq+=ops_set[stemSeq[-j]]
                    else:
                        altSeq+=random.choice(nt_set)
                hp_seq=altSeq+apicalSeq+stemSeq

            else:
                apicalSeq=fullSeq[len(fullSeq)-apicalSize:]
                stemSeq=fullSeq[:len(fullSeq)-apicalSize]
                for j in range(1,len(stemSeq)+1):
                    check = random.randint(0,100)
                    if check <= paired_percent or j > len(stemSeq)-3:
                        altSeq+=ops_set[stemSeq[-j]]
                    else:
                        altSeq+=random.choice(nt_set)
                hp_seq=stemSeq+apicalSeq+altSeq

            tick=""
            if args.Repeats:
                tick+=str(repWindow)+"-"
                foundRepeat=False
                repList=[]
                for i in range(0, len(hp_seq)):
                    subRep=hp_seq[i:i+repWindow]
                    if len(subRep) >= repWindow:
                        p=[match.start() for match in re.finditer(subRep, hp_seq)]
                        if len(p) > 1:
                            foundRepeat=True
                            if subRep not in repList:
                                tick+=subRep+":"
                                repList.append(subRep)
                            for t in p:
                                if t not in repList:
                                    repList.append(t)
                                    tick+=str(t+1)+":"
                            tick=tick[0:len(tick)-1]
                            tick+="_"
                if len(tick) > 2:
                    tick=tick[0:len(tick)-1]

            if os.path.isfile(outFile):
                data = ""
            else:
                data = "Name,Complementarity,ApicalSize,ApicalSeq,Bulge,BulgeSize,BulgePosition,BulgeSeq,StemLength,Sequence,Structure,Length,bp,GC,dG,dG_Length,PE,APE,GC_pairs,AU_pairs,GU_pairs,GC_pair_percent,AU_pair_percent,GU_pair_percent,apicals,left_bulges,right_bulges,repeats\n"

            fc = RNA.fold_compound(hp_seq)
            fc.pf()
            
            pairs = get_pairs(fc.mfe()[0], hp_seq)
            bulge_counts = bulge_count(fc.mfe()[0], hp_seq)

            data+=i[1:]+"_"+aType+"_"+str(iter+1)+"," #Adds a general name
            data+=str(paired_percent)+"," #Adds complementarity score
            data+=str(apicalSize)+"," # Adds the size of the apical loop
            data+=apicalSeq+"," # Adds the apical sequence
            data+="None," #Adds the bulge type
            data+="0," #Adds Bulge size
            data+="None," #Add Bugle position
            data+="None," #Adds bulge sequence
            data+=str(len(stemSeq))+"," #Adds the length of hairpin stem
            data+=hp_seq+"," #Adds the primary sequence
            data+=fc.mfe()[0]+"," # Adds the structure sequence in dot-bracket
            data+=str(len(hp_seq))+"," # Adds the total length of the sequence
            data+=str(fc.mfe()[0].count("("))+"," # Adds the number of basepaired bases
            data+=str(round((hp_seq.count("G")+hp_seq.count("C"))/len(hp_seq),2))+"," # Adds the GC content
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
            if args.Details:
                if len(bulge_counts["apical_local"]) > 0:
                    data+=str(bulge_counts["apical_local"][0:len(bulge_counts["apical_local"])-1])+"," #Adds the number of apical loops
                else:
                    data+="None,"
                if len(bulge_counts["lB_local"]) > 0:
                    data+=str(bulge_counts["lB_local"][0:len(bulge_counts["lB_local"])-1])+"," #Adds the number of 5 prime bulges
                else:
                    data+="None,"
                if len(bulge_counts["rB_local"]) > 0:
                    data+=str(bulge_counts["rB_local"][0:len(bulge_counts["rB_local"])-1])+"," #Adds the number of 3 prime bulges
                else:
                    data+="None,"
            else:
                data+=str(bulge_counts["apicals"])+"," #Adds the number of apical loops
                data+=str(bulge_counts["l_bulges"])+"," #Adds the number of 5 prime bulges
                data+=str(bulge_counts["r_bulges"])+"," #Adds the number of 3 prime bulges
            if len(tick) > 2:
                if args.Details:
                    data+=tick+"\n" #Tags if sequence contains a repeat
                else:
                    data+=str(repWindow)+"-"+str(foundRepeat)+"\n" #Tags if sequence contains a repeat
            elif args.Repeats:
                data+=str(repWindow)+"-"+str(foundRepeat)+"\n" #Tags if sequence contains a repeat
            else:
                data+="NA\n" #Tags if sequence contains a repeat

            with open(outFile, 'a') as f:
                f.write(data)

    print("\nHairpins generated and stored in "+outFile)