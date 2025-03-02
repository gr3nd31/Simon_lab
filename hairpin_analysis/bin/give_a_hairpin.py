import argparse
import RNA
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file. If blank, a random hairpin is generated") #
parser.add_argument("-l", "--length", help="Maximum length sequence to be folded")
parser.add_argument("-R", "--Repeats", help="Flag repeats. Give an integer length to check if sequence is repeated.")
parser.add_argument("-D", "--Details", help="Report sequences of apical loops and bulges rather than number.",  action='store_true')
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()
sliceIt = False

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

# Names the output file
if args.out:
    outFile = args.out
else:
    outFile = "data.csv"
# Pulls of the given sequence
if args.sequence:
    seqs = read_fasta(args.sequence)
    try:
        seqs = read_fasta(args.sequence)
    except:
        print("Unable to parse given fasta file.")
        seqs = "nope"
else:
    print("Fasta file not given. Use '-s' flag and path to file.")
    seqs = "nope"
# Triggers the string to be split in n-sized chunks, best used for very large RNA molecules
if args.length:
    sliceIt=True
    slicer=int(args.length)

if seqs != "nope":
    for iter in seqs:
        hp_seq=seqs[iter].upper().replace(" ", "")
        hp_seq=hp_seq.replace("T", "U")
        if sliceIt:
            slices = [hp_seq[a:a+slicer] for a in range(0, len(hp_seq), slicer)]
#            if len(slices[len(slices)-1]) < slicer:
#                slices[len(slices)-1] = hp_seq[len(hp_seq)-slicer:len(hp_seq)]
        else:
            slices = [hp_seq]
        
        for subSeq in slices:

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
                            for t in p:
                                if t not in repList:
                                    repList.append(t)
                                    tick+=str(t+1)+":"
                            tick=tick[0:len(tick)-1]
                            tick+="_"
                if len(tick) > 2:
                    tick=tick[0:len(tick)-1]

            fc = RNA.fold_compound(subSeq)
            fc.pf()

            pairs = get_pairs(fc.mfe()[0], subSeq)
            bulge_counts = bulge_count(fc.mfe()[0], subSeq)

            structure = fc.mfe()[0].replace(",", ".")
            paired_percent = 1-round(structure.count(".")/len(structure),2)
            #print(paired_percent)
            ap_length = 0
            ap_seq = ""
            codpiece = "None"
            baseBulge = ""
            where_the_b = "0.5"
            hp_length = structure.count("(")

            if os.path.isfile(outFile):
                data = ""
            else:
                data = "Name,Complementarity,ApicalSize,ApicalSeq,Bulge,BulgeSize,BulgePosition,BulgeSeq,StemLength,Sequence,Structure,Length,bp,GC,dG,dG_Length,PE,APE,GC_pairs,AU_pairs,GU_pairs,GC_pair_percent,AU_pair_percent,GU_pair_percent,apicals,left_bulges,right_bulges,repeats\n"
            
            data+=iter.replace("<", "").replace(",", "")+"," #Adds a general name and removes any devilish commas
            data+=str(paired_percent)+"," #Adds complementarity score
            data+=str(ap_length)+"," # Adds the size of the apical loop
            data+=ap_seq+"," # Adds the apical sequence
            data+=codpiece+"," #Adds the bulge type
            data+=str(len(baseBulge))+"," #Adds Bulge size
            data+=where_the_b+"," #Add Bugle position
            data+=baseBulge+"," #Adds bulge sequence
            data+=str(hp_length)+"," #Adds the length of hairpin stem
            data+=subSeq+"," #Adds the primary sequence
            data+='"'+fc.mfe()[0]+'",' # Adds the structure sequence in dot-bracket
            data+=str(len(subSeq))+"," # Adds the total length of the sequence
            data+=str(fc.pf()[0].count("("))+"," # Adds the number of basepaired bases
            data+=str(round((subSeq.count("G")+subSeq.count("C"))/len(subSeq),2))+"," # Adds the GC content
            data+=str(round(fc.mfe()[1],3))+"," #Adds the MFE
            data+=str(round(fc.mfe()[1]/len(subSeq),2))+"," # Adds the MFE/length ratio
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
                data+=str(bulge_counts["apical_local"])+"," #Adds the number of apical loops
                data+=str(bulge_counts["lB_local"])+"," #Adds the number of 5 prime bulges
                data+=str(bulge_counts["rB_local"])+"," #Adds the number of 3 prime bulges
            else:
                data+=str(bulge_counts["apicals"])+"," #Adds the number of apical loops
                data+=str(bulge_counts["l_bulges"])+"," #Adds the number of 5 prime bulges
                data+=str(bulge_counts["r_bulges"])+"," #Adds the number of 3 prime bulges
            if len(tick) > 2:
                data+=tick+"\n" #Tags if sequence contains a repeat
            elif args.Repeats:
                data+=str(repWindow)+"-"+str(foundRepeat)+"\n" #Tags if sequence contains a repeat
            else:
                data+="NA\n" #Tags if sequence contains a repeat

            with open(outFile, 'a') as f:
                f.write(data)