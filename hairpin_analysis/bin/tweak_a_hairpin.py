import argparse
import RNA
import random
import os
import re
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--entropy", help = "Requested final delta G/length ratio (default = -0.43)", default=-0.43)
parser.add_argument("-n", "--number", help = "Number of iterations per sequence.", default=1)
parser.add_argument("-A", "--ApicalLength", help="Length of the apical sequence.", default=5)
parser.add_argument("-R", "--Repeats", help="Flag repeats. Give an integer length to check if sequence is repeated.", default=6)
parser.add_argument("-s", "--sequence", help="Path to the fasta to be tweaked. If sequence doesn not form a hairpin, make_a_hairpin.py is called to generate a hairpin from the given sequence.")
parser.add_argument("-a", "--max_APE", help="Maximum APE value (Default is 0.4)", default=0.4)
parser.add_argument("-p", "--max_paired_length", help="Maximum length of perfectly paired RNA.", default=21)
parser.add_argument("-c", "--coding", help="If hairpin has coding, gives the position at which a codon begins and attempts to only manipulate nucleotides that maintain coding.")
parser.add_argument("-o", "--out", help="Name of output file (Default is 'data.csv')", default="data.csv")
parser.add_argument("-u", "--uncertainty", help="Percent (0-1) of target delta G/Length sufficient for tweaking.", default=0.1)
parser.add_argument("-F", "--Force_hairpin", help="If flagged, sequences always used as the 5' side of a generated hairpin.",  action='store_true')
parser.add_argument("-C", "--ConserveSequence", help="If flagged, the input sequence cannot be mutated.")
parser.add_argument("-S", "--Subset", help="List of positions that should be changed. Example format: 12,15,20-30,60-70")
parser.parse_args()
args = parser.parse_args()
runit = True

# Set to randomly pull from
nt_set = ["A", "U", "G", "C"]
# Dictionary to find complements
ops_set = {"A":"U", "U":"A", "C":"G", "G":"C"}
# Set of positions that cannot be changed
unFunSet = []

# Codon usage table
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

# Get positions of paired partner
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

# Function to count apical loops and bulges
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

# Fasta reading function
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

    if len(sequences) > 1:
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
    elif len(sequences) == 1:
        print("FASTA file not given. Assuming first line is a sequence and using placeholder name of 'sequence'")
        seq_dic_2 = {"sequence": sequences[0]}

    return seq_dic_2

# Function to fold RNA and return variables
def fold_and_return(sequence):
    folded = RNA.fold_compound(sequence)
    folded.pf()
    structure = folded.mfe()[0]
    current_dg = folded.mfe()[1]
    pes = list(folded.positional_entropy())[1:]
    ape = sum(pes)/len(pes)
    return structure, current_dg, pes, ape

# Splits structure into paired partners
def split_and_knit(sequence, structure):
    box = {}
    paired=[]
    unpairs=[]
    bCount = 1
    real_start=1
    real_end=len(sequence)
    max_pair = 0
    pair_length = 0
    while structure.find("(") >= 0:
        start = structure.find("(")
        end = structure.rfind(")")
        if start > 0 or end < len(sequence)-1:
            box[bCount] = ["U", real_start, real_end, sequence[:start],sequence[end+1:]]
            unpairs.append(bCount)
            bCount+=1
            real_start+=len(sequence[:start])
            real_end-=len(sequence[end+1:])
            if pair_length > max_pair:
                max_pair=pair_length
            pair_length = 1
        else:
            pair_length+=1
        box[bCount] = ["P", real_start, real_end, sequence[start],sequence[end]]
        real_start+=1
        real_end-=1
        paired.append(bCount)
        bCount+=1
        structure = structure[start+1:end]
        sequence = sequence[start+1:end]

    box[bCount] = ["A", real_start, real_end, sequence, ""]
    rCounter = 0
    apical = sequence
    unpairs.append(bCount)
    return box, apical, rCounter, max_pair

# Function to find repeat sequences
def repeat_finder(sequence, repWindow):
    repList={}
    for r in range(0, len(sequence)):
        subRep=sequence[r:r+repWindow]
        if len(subRep) >= repWindow:
            p=[match.start() for match in re.finditer(subRep, sequence)]
            if len(p) > 1:
                if subRep not in repList.keys():
                    repList[subRep] = p
    return repList

# Makes dataframe of the positions, paired partners, and possible changes/codons
def new_split(sequence, box, cons, codes):
    fox = {"position":[],"base":[],"paired":[],"partner":[], "set":[],"fungible":[], "possibles":[], "codonCons":[]}
    counter=1
    for i in sequence:
        fox["position"].append(counter)
        fox["base"].append(i)
        fox["paired"].append("U")
        fox["partner"].append("test")
        fox["set"].append(0)
        if counter-1 in cons.keys():
            #print(counter-1)
            fox["fungible"].append("Yes")
        else:
            fox["fungible"].append("No")
        fox["possibles"].append("Hold")
        fox["codonCons"].append("Hold")
        counter+=1
    box_df = pd.DataFrame(fox)
    for i in cons:
        ps=""
        for j in cons[i]:
            ps+=j+"_"
        ps = ps[:len(ps)-1]
        box_df.loc[box_df['position'] == i+1, 'possibles'] = ps

    if len(codes) > 0:
        for i in codes:
            a_list=""
            amino=i[2]
            for a in amino:
                a_list+="_"+a
            starter=i[1][0]
            box_df.loc[box_df['position'] == starter+1, 'codonCons'] = "1"+a_list
            box_df.loc[box_df['position'] == starter+2, 'codonCons'] = "2"+a_list
            box_df.loc[box_df['position'] == starter+3, 'codonCons'] = "3"+a_list

    for i in box:
        interim = box[i]
        if interim[0] == "P":
            box_df.loc[box_df['position'] == interim[1], 'partner'] = i
            box_df.loc[box_df['position'] == interim[2], 'partner'] = i
            if (interim[3] == "G" or interim[4] == "G") and (interim[3] == "C" or interim[4] == "C"):
                box_df.loc[box_df['position'] == interim[1], 'paired'] = "GC"
                box_df.loc[box_df['position'] == interim[2], 'paired'] = "GC"
                box_df.loc[box_df['position'] == interim[1], 'set'] = 1
                box_df.loc[box_df['position'] == interim[2], 'set'] = 2
            elif (interim[3] == "A" or interim[4] == "A") and (interim[3] == "U" or interim[4] == "U"):
                box_df.loc[box_df['position'] == interim[1], 'paired'] = "AU"
                box_df.loc[box_df['position'] == interim[2], 'paired'] = "AU"
                box_df.loc[box_df['position'] == interim[1], 'set'] = 1
                box_df.loc[box_df['position'] == interim[2], 'set'] = 2
            elif (interim[3] == "G" or interim[4] == "G") and (interim[3] == "U" or interim[4] == "U"):
                box_df.loc[box_df['position'] == interim[1], 'paired'] = "GU"
                box_df.loc[box_df['position'] == interim[2], 'paired'] = "GU"
                box_df.loc[box_df['position'] == interim[1], 'set'] = 1
                box_df.loc[box_df['position'] == interim[2], 'set'] = 2
        else:
            if len(interim[3]) > 0:
                for j in range(0,len(interim[3])):
                    box_df.loc[box_df['position'] == interim[1]+j, 'partner'] = i
                    box_df.loc[box_df['position'] == interim[1]+j, 'set'] = 1
            if len(interim[4])> 0:
                for j in range(0,len(interim[4])):
                    box_df.loc[box_df['position'] == interim[2]-j, 'partner'] = i
                    box_df.loc[box_df['position'] == interim[2]-j, 'set'] = 2
            
    return box_df

# Defines how close the delta G should be from the desired parameter
if args.uncertainty:
    try:
        threshold=float(args.uncertainty)
    except:
        print("Unable to read given uncertainty. Aborting")
        runit=False

# Defines the number of iterations to return
if args.number:
    try:
        replNum=int(args.number)
    except:
        print("Unable to parse given integer. Defaulting to 1")
        replNum=1

# Option to check is a sequence is repeated in the hairpin
if args.Repeats:
    try:
        repWindow=int(args.Repeats)
        foundRepeat=False
    except:
        print("Unable to parse given integer. Defaulting to 6.")
        repWindow=6
        foundRepeat=False
else:
    repWindow=0
    foundRepeat=False

# Define the target delta G
if args.entropy:
    try:
        target_dG = float(args.entropy)
        if target_dG > 0 or target_dG < -1:
            print("Given parameters are impossible for a hairpin (>0 or less than < -1). Aborting")
            runit=False
    except:
        print("Unable to read given entropy value. Aborting.")
        runit=False

# Pulls the fasta file
if args.sequence:
    try:
        seqs = read_fasta(args.sequence)
    except:
        print("Unable to parse given fasta file.")
        runit=False
else:
    print("Fasta file not given. Use '-s' flag and path to file.")
    runit=False

#Define the maximum APE
if args.max_APE:
    try:
        max_ape = float(args.max_APE)
        if max_ape < 0:
            print("Given APE is not possible (below 0). Aborting")
            runit=False
    except:
        print("Unable to read given APE value. Aborting")
        runit=False

# Define the size of the apical loop OR the sequence to be used in the apical loop
if args.ApicalLength:
    ap_seq=""
    try:
        ap_length=int(args.ApicalLength)
    except:
        tryIt=input("Apical input is detected as a string. Would you like to use this as the apical sequence (y/N)? ")
        if tryIt.upper() == "Y":
            ap_seq = args.ApicalLength
            ap_seq = ap_seq.upper().replace('T', 'U')
            ap_length=len(ap_seq)
        else:
            print("Disregarding the string input and using default (5) apical length.")
            ap_length=5

# Define value of maximum number to consecutive base pairs.
if args.max_paired_length:
    try:
        max_paired = int(args.max_paired_length)
    except:
        print("Unable to read given max paired length. Aborting")
        runit = False

# Flag to tell script to conserve coding beginning at a certain position
if args.coding:
    try:
        codon_start = int(args.coding)
        if codon_start < 1:
            print("Given codon cannot be used. Aborting.")
            runit=False
        else:
            print("\nUsing standard codon table to identify which positions can be changed starting at position "+args.coding+".")
    except:
        print("Unable to read given codon start. Aborting.")
        runit=False

# Flag to tell script to conserve the input sequence
if args.ConserveSequence:
    try:
        unfunPos=args.ConserveSequence.split(",")
        for j in unfunPos:
            thing=j.split("-")
            if len(thing) == 1:
                unFunSet.append(int(thing[0])-1)
            elif len(thing) == 2:
                for i in list(range(int(thing[0]), int(thing[1])+1)):
                    unFunSet.append(i-1)
    except:
        print("Unable to parse the given conserved positions. Aborting.")
        runit=False

# Reads list of positions that CAN be changed
if args.Subset:
    try:
        aList =[]
        ks = args.Subset.split(",")
        for i in ks:
            if "-" in i:
                z=i.split("-")
                for j in range(int(z[0])-1, int(z[1])):
                    aList.append(j)
            else:
                aList.append(int(i))
    except:
        print("Unable to read input list. Try again using the correct format. Aborting")
        runit=False

# Assuming all the inputs are correct, the script runw
if runit:
    # Iterates through the fasta sequences
    for i in seqs:
        # Repeats for the number of requested iterations
        for n in range(0,replNum):
            # Defines an initial blank apical sequence
            apical = ""
            # Tells user which sequence and iteration is being run
            print("\nWorking on sequence: "+i[1:]+" iteration "+str(n+1))
            # Converts sequence to RNA and makes uppercase
            working_seq = seqs[i].upper().replace('T', 'U')
            # Gets initial structure and thermodynamic values to see if the sequence is *already* a hairpin
            structure, current_dg, pes, ape = fold_and_return(working_seq)
            # Count apicals to check its a hairpin
            bulge_counts = bulge_count(structure, working_seq)
            # if there is not a single apical loop or if the user has requested to make a new hairpin from the sequence, a new hairpin is generated around the given sequence
            if (bulge_counts["apicals"] != 1 or args.Force_hairpin) and not args.coding:
                # Defines a range of acceptable delta G/Length values
                target_G = ((len(working_seq)*2)+5)*target_dG
                max_allowed = target_G+(target_G*threshold)
                min_allowed = target_G-(target_G*threshold)
                # Initiates a variabel to determine when all parameter threshold are met
                allGood = False
                # Alerts user that a new hairpin is being generated
                if bulge_counts["apicals"] != 1 and not args.Force_hairpin:
                    print("Given sequence is not predicted to form a hairpin. Using sequence to form a hairpin.")
                else:
                    print("Duplicating sequence to form a hairpin.")
                # Estimates a percent complementarity to achieve a given deltaG/Length based on the linear regression of deltaG/Length and percent complementarity
                basecompl = round((target_dG-0.248)/-1.044580, 2)
                print("Based on the target dG/Length ratio, the hairpin should have a "+str(basecompl*100)+"% complementarity.")
                # Initiates counter variable so the user can be alerted if a given sequence should be stopped
                counter = 0
                # while loop until all parameters are met
                while not allGood:
                    # Refreshes the working sequence
                    working_seq = seqs[i].upper().replace('T', 'U')
                    # Generates a new blank apical loop
                    apical=""
                    # Generates a new complement sequence
                    complement=""
                    # If a specific apical sequence has not been requested, it is randomly generated
                    if ap_seq == "":
                        for j in range(0, ap_length):
                            apical+=random.choice(nt_set)
                    else:
                        apical+=ap_seq

                    # A complementary sequence is generated that randomly generate mismatches based on the calculated complementarity
                    for j in working_seq[::-1]:
                        # Randomize things based on estimated probability
                        check = random.randint(0,100)
                        if check <= basecompl*100:
                            complement+=ops_set[j]
                        else:
                            possibility=random.choice(nt_set)
                            while possibility == j:
                                possibility=random.choice(nt_set)
                            complement+=possibility
                    # assembles the hairpin sequence
                    working_seq+=apical+complement
                    # Folds the generated hairpin to get structure and thermodynamic values
                    structure, current_dg, pes, ape = fold_and_return(working_seq)
                    # Records apical loop and bulge information
                    bulge_counts = bulge_count(structure, working_seq)
                    # If the thermodynamic values are within the defined limits, shows a strongly paired base, and only has a single apical loop, the sequence is completed and the loop is broken
                    if bulge_counts["apicals"] == 1 and current_dg >= max_allowed and current_dg <= min_allowed and structure.startswith("((((") and ape < max_ape:
                        allGood = True
                    # Increase counter
                    counter+=1
                    # If iteration has occured at increments of 500, the user is given the chance to break the loop manually
                    if counter%500 == 0:
                        cont = input("Attempted "+str(counter)+" iterations to match delta G. Continue (y/n)? ")
                        if cont.lower() == "n":
                            allGood=True

            # Thermodynamic parameters are reported to the user
            target_G = len(working_seq)*target_dG
            max_allowed = target_G+(target_G*threshold)
            min_allowed = target_G-(target_G*threshold)
            print("Target delta G: "+str(round(target_G,2))+" ("+str(round(min_allowed, 2))+" - "+str(round(max_allowed, 2))+", ratio:"+str(target_dG)+" delta G/length)")

            # Information about the pairing of each base is obtained
            box, apical, rCounter, maxPairs = split_and_knit(working_seq, structure)
            subSeq = working_seq

            # Identify which positions can be changed
            fungible_positions = {}
            condonConsiderations=[]
            # If coding must be maintained, the potential positions that can be altered are flagged
            if args.coding:
                # iterate through codons
                for p in range(codon_start-1,len(working_seq),3):
                    codon=working_seq[p:p+3]
                    if len(codon) == 3:
                        # Get the possible codons for that amino acid
                        keys = [key for key, val in codon_table.items() if val == codon_table[codon]]
                        condonConsiderations.append([codon_table[codon],[p, p+1, p+2], keys])
                        # If the number of possible codons is greater than one, the possibilities are pulled
                        if len(keys) > 1:
                            # Possible bases of each position are pulled
                            for c in range(0,3):
                                for k in keys:
                                    if k[c] != codon[c]:
                                        if p+c+1 in fungible_positions.keys():
                                            fungible_positions[p+c+1].append(k[c])
                                        else:
                                            fungible_positions[p+c] = [k[c]]
                                if p+c in fungible_positions.keys():
                                    fungible_positions[p+c].append(codon[c])
                                    fungible_positions[p+c] = list(set(fungible_positions[p+c]))
            
            # If given sequence must be conserved, the input positions are flagged to not be changed
            elif args.ConserveSequence:
                funNums = list(range(1,len(working_seq)+1))
                for c in funNums:
                    if c not in unFunSet:
                        fungible_positions[c]=nt_set
            # If specific positions are flagged as changeable, those positions are added the list of positions to be altered
            elif args.Subset:
                funNums = list(range(0,len(working_seq)))
                for c in funNums:
                    if c in aList:
                        fungible_positions[c]=nt_set
            # otherwise, all lists are considered fungible
            else:
                for c in list(range(0, len(working_seq))):
                    fungible_positions[c]=nt_set

            # Tells user the initial thermodynamic parameters
            print("Initial delta G is "+str(round(current_dg,2))+" (dG/Length: "+str(round(current_dg/len(subSeq),2))+" kilocalories/mole/nt.)")
            # If the desired values are not yet met, the user is altered
            if current_dg <= max_allowed or current_dg >= min_allowed:
                print("Initial hairpin requires delta G tweaking.")
                tweak_dG = True
            else:
                #print("Initial hairpin fits delta G, moving to PE.")
                tweak_dG = False
            
            df = new_split(sequence=subSeq, box=box, cons=fungible_positions, codes=condonConsiderations)
#----------------------------------------------------------------------------- Tweak delta G
            # now we change pairing to get delta G right
            if tweak_dG:
                # Counter initiated/reset
                counter = 0
                ticket = 0
                initTry = 100
                # while loop until deltaG is allowed or the structure is no longer a hairpin
                while current_dg <= max_allowed or current_dg >= min_allowed or bulge_counts["apicals"] != 1:
                    # Variable to track when changes are made
                    changed=False
                    # Dictionary of pairs/unpaired positions are generated
                    # This allows the script to make or break bonds to tweak teh delta G
                    gc={}
                    au={}
                    gu={}
                    unpaired={}
                    for b in box:
                        if box[b][0] == 'P':
                            if (box[b][3] == "A" or box[b][4] == "A") and (box[b][3] == "U" or box[b][4] == "U"):
                                au[b]=box[b]
                            elif (box[b][3] == "G" or box[b][4] == "G") and (box[b][3] == "C" or box[b][4] == "C"):
                                gc[b]=box[b]
                            else:
                                gu[b]=box[b]
                        else:
                            if len(box[b][3]) > 0 and len(box[b][4]) > 0:
                                unpaired[b]=box[b]

                    # If the deltaG is too low, basepaired positions are selected
                    if current_dg <= max_allowed:
                        trix = df.loc[(df.fungible == "Yes") & (df.paired != "U")]
                        if len(trix) > 0:
                            trix = trix[trix.position == random.choice(list(trix.position))]
                    # If the deltaG is too high, unpaired positions are selected
                    elif current_dg >= min_allowed:
                        trix = df.loc[(df.fungible == "Yes") & (df.paired == "U")]

                    while len(trix) > 0:
                        # Get random unpaired position that share a box
                        # If coding, the codon is changed
                        if args.coding:
                            trix = trix[trix.position == random.choice(list(trix.position))]
                            subSeq=""
                            for x in list(trix.position):
                                possibles = random.choice(trix.loc[trix.position == x, 'codonCons'].item().split("_")[1:])
                                offset = int(trix.loc[trix.position == x, 'codonCons'].item().split("_")[0])
                                for b in range(0,3):
                                    gym=(x-offset)+b+1
                                    df.loc[df.position == gym, 'base'] = possibles[b]
                                trix = trix[trix.position != x]
                                changed=True
                        # Otherwise, the base is changed randomly
                        else:
                            subSeq=""
                            for x in list(trix.position):
                                possibles = trix.loc[trix.position == x, 'possibles'].item().split("_")
                                df.loc[df.position == x, 'base'] = random.choice(possibles)
                                trix = trix[trix.position != x]
                            changed=True
                        # The sequence is put back together from the dataframe
                        for b in range(1,max(df.position)+1):
                            subSeq+=df.loc[df.position == b, 'base'].item()

                    # New sequence folded and pairing partners re-evaluated
                    structure, current_dg, pes, ape = fold_and_return(subSeq)
                    box, apical, rCounter, maxPairs = split_and_knit(subSeq, structure)
                    bulge_counts = bulge_count(structure, subSeq)
                    df = new_split(sequence=subSeq, box=box, cons=fungible_positions, codes=condonConsiderations)
                    counter+=1

                    # User is allowed to break the loop every 100 iterations
                    if counter >= initTry and bulge_counts["apicals"] == 1:
                        ticket+=counter
                        keepIt = input("Have tried "+str(ticket)+" permutations. The current delta G is at "+str(round(current_dg, 2))+". Should the function keep trying (y/N)? ")
                        if keepIt.upper() == "Y":
                            counter = 0
                        else:
                            break
                    # If the mutations create a branched structure, the user can 'reset' the sequence 
                    elif counter >= initTry and bulge_counts["apicals"] > 1:
                        ticket+=counter
                        keepIt = input("Have tried "+str(ticket)+" permutations. The current delta G is at "+str(round(current_dg, 2))+" but multiple apical loops "+str(bulge_counts["apicals"])+" detected. Should the sequence be reset to the original (y/N)? ")
                        if keepIt.upper() == "Y":
                            subSeq = working_seq
                            structure, current_dg, pes, ape = fold_and_return(subSeq)
                            box, apical, rCounter, maxPairs = split_and_knit(subSeq, structure)
                            bulge_counts = bulge_count(structure, subSeq)
                            df = new_split(sequence=subSeq, box=box, cons=fungible_positions, codes=condonConsiderations)
                        else:
                            keepIt = input("Would you like to continue with the current sequence despite having multiple apical loops? (y/N) ")
                            if keepIt.upper() != "Y":
                                break
                        counter = 0
                    # if thermodynamics are achieved, the loop is broken
                    else:
                        if current_dg >= max_allowed and current_dg <= min_allowed and bulge_counts["apicals"] == 1:
                            #print("Appropriate delta G has been achieved. Checking PE values.")
                            break
                print("Appropriate delta G has been achieved ("+str(round(current_dg, 2))+"). Checking PE values.\n")
            else:
                print("Delta G already appropriate. Checking PE values.\n")
#----------------------------------------------------------------------------- Tweak PE
            # If the APE is too high, the highest PE bases are mutated
            if ape > max_ape:
                peGO = input("Average postitional entropy is too high ("+str(round(ape,2))+"). Would you like to try decreasing PE (CAUTION: CHANGING HIGH PE BASES MAY CHANGE THE DELTA G)? (y/N) ")
                if peGO.upper() == "Y":
                    if args.coding:
                        print("Sorry, altering PE in a coding sequence is too complex. Try adjusting delta G.")
                    else:
                        counter=0
                        ticket=0
                        while ape > max_ape:
                            # The highest PE base is identified
                            max_pos = pes.index(max(pes))
                            # if the highest PE base cannot be changed, it is ignored
                            if max_pos+1 in unFunSet:
                                pes.remove(max(pes))
                            else:
                                possibles = df.loc[df.position == max_pos+1, 'possibles'].item().split("_")
                                pairing = df.loc[df.position == max_pos+1, 'paired'].item()
                                # If the highest PE base is unpaired, it is randomly mutated
                                if pairing == "U":
                                    df.loc[df.position == max_pos+1, 'base'] = random.choice(possibles)
                                # If the highest PE base is paired, the base pairing is mutated
                                else:
                                    partner = df.loc[df.position == max_pos+1, 'partner'].item()
                                    if pairing == "GC":
                                        pickIt = random.choice(["AU", "UA"])
                                        df.loc[df.position == max_pos+1, 'base'] = pickIt[0]
                                        df.loc[df.position == partner, 'base'] = pickIt[1]
                                    elif pairing == "AU":
                                        pickIt = random.choice(["GC", "CG"])
                                        df.loc[df.position == max_pos+1, 'base'] = "C"
                                        df.loc[df.position == max_pos+1, 'base'] = pickIt[0]
                                        df.loc[df.position == partner, 'base'] = pickIt[1]
                                    else:
                                        if subSeq[max_pos == "U"]:
                                            df.loc[df.position == max_pos+1, 'base'] = "C"
                                        else:
                                            df.loc[df.position == max_pos+1, 'base'] = "A"
                                
                                # Sequence is reconstructed
                                subSeq=""
                                for b in range(1,max(df.position)+1):
                                    subSeq+=df.loc[df.position == b, 'base'].item()
                                # New sequence is folded and the APE is evaluated
                                structure, current_dg, pes, ape = fold_and_return(subSeq)
                                box, apical, rCounter, maxPairs = split_and_knit(subSeq, structure)
                                bulge_counts = bulge_count(structure, subSeq)
                                df = new_split(sequence=subSeq, box=box, cons=fungible_positions, codes=condonConsiderations)
                                # Counter to allow the user to escape the loop
                                if counter > 100:
                                    ticket+=counter
                                    keepIt = input("Have tried "+str(ticket)+" permutations. The current APE is at "+str(round(ape, 2))+". Should the function keep trying (y/N)? ")
                                    if keepIt.upper() == "Y":
                                        counter = 0
                                    else:
                                        break
                                counter+=1
                        print("Exiting PE adjustment. Current APE is "+str(round(ape,2))+".")
                        if current_dg >= max_allowed and current_dg <= min_allowed and bulge_counts["apicals"] == 1:
                            print("Delta G of the hairpin is still within parameters ("+str(round(current_dg, 2))+"). All parameters met.\n")
                        else:
                            print("Delta G of the hairpin ("+str(round(current_dg, 2))+") is now outside of acceptable range. Try another iteration.\n")
            else:
                print("APE ("+str(round(ape,2))+") already below target threshold ("+str(max_ape)+"). Thermodynamic parameters met.\n")
#----------------------------------------------------------------------------- Reports if the pairing is too great
            if maxPairs > args.max_paired_length:
                print("Detected a perfectly paired length ("+str(maxPairs)+") greater than allowed. Consider decreasing delta G parameter.\n")
            else:
                print("Perfectly paired length ("+str(maxPairs)+") is below the given threshold ("+str(args.max_paired_length)+").\n")

#----------------------------------------------------------------------------- Tweak repeats
            repeaters = repeat_finder(subSeq, repWindow)
            tick = ""
            if args.Repeats:
                tick+=str(repWindow)+"-"
                if len(repeaters) > 0:
                    print("Found "+str(len(repeaters))+" repeats found of at least "+str(args.Repeats)+" bases long.")
                    for r in repeaters:
                        tick+=r+":"
                        for s in repeaters[r]:
                            tick+=str(s)+":"
                        tick = tick[0:len(tick)-1]
                        tick+="_"
                    tick = tick[0:len(tick)-1]
                else:
                    print("No repeats found of "+str(args.Repeats)+" bases long.")
                    foundRepeat=False
    
#-----------------------------------------------------------------------------
            print("Writing hairpin data to file")
            # Final sequence is folded and analyzed
            fc = RNA.fold_compound(subSeq)
            fc.pf()
            pairs = get_pairs(fc.mfe()[0], subSeq)
            bulge_counts = bulge_count(fc.mfe()[0], subSeq)

            structure = fc.mfe()[0].replace(",", ".")
            paired_percent = 1-round(structure.count(".")/len(structure),2)
            codpiece = "None"
            baseBulge = ""
            where_the_b = "0.5"
            hp_length = structure.count("(")

            if os.path.isfile(args.out):
                data = ""
            else:
                data = "Name,Complementarity,ApicalSize,ApicalSeq,Bulge,BulgeSize,BulgePosition,BulgeSeq,StemLength,Sequence,Structure,Length,bp,GC,dG,dG_Length,PE,APE,GC_pairs,AU_pairs,GU_pairs,GC_pair_percent,AU_pair_percent,GU_pair_percent,apicals,left_bulges,right_bulges,repeats\n"
                
            data+=i.replace("<", "").replace(",", "")+"_"+str(n)+"," #Adds a general name and removes any devilish commas
            data+=str(paired_percent)+"," #Adds complementarity score
            data+=str(len(apical))+"," # Adds the size of the apical loop
            data+=apical+"," # Adds the apical sequence
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

            if len(tick) > 2:
                data+=tick+"\n" #Tags if sequence contains a repeat
            else:
                data+="NA\n" #Tags if sequence contains a repeat

            with open(args.out, 'a') as f:
                f.write(data)
    print("Hairpin design complete. See "+args.out+" for sequence details.")
        