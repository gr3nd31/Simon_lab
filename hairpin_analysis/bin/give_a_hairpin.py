import argparse
import RNA
import os

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file. If blank, a random hairpin is generated") #
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()

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

# Names the output file
if args.out:
    outFile = args.out
else:
    outFile = "data.csv"

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

seqs = read_fasta(args.sequence)
for iter in seqs:
    #print(iter)
    hp_seq=seqs[iter].upper()
    #print(hp_seq)
    fc = RNA.fold_compound(hp_seq)
    fc.pf()

    pairs = get_pairs(fc.mfe()[0], hp_seq)

    structure = fc.mfe()[0].replace(",", ".")
    paired_percent = round(structure.count(".")/len(structure),2)
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
        data = "Name,Complementarity,ApicalSize,ApicalSeq,Bulge,BulgeSize,BulgePosition,BulgeSeq,StemLength,Sequence,Structure,Length,bp,GC,dG,dG_Length,PE,APE,GC_pairs,AU_pairs,GU_pairs,GC_pair_percent,AU_pair_percent,GU_pair_percent\n"

    data+=iter.replace("<", "")+"," #Adds a general name
    data+=str(paired_percent)+"," #Adds complementarity score
    data+=str(ap_length)+"," # Adds the size of the apical loop
    data+=ap_seq+"," # Adds the apical sequence
    data+=codpiece+"," #Adds the bulge type
    data+=str(len(baseBulge))+"," #Adds Bulge size
    data+=where_the_b+"," #Add Bugle position
    data+=baseBulge+"," #Adds bulge sequence
    data+=str(hp_length)+"," #Adds the length of hairpin stem
    data+=hp_seq+"," #Adds the primary sequence
    data+='"'+fc.mfe()[0]+'",' # Adds the structure sequence in dot-bracket
    data+=str(len(hp_seq))+"," # Adds the total length of the sequence
    data+=str(fc.pf()[0].count("("))+"," # Adds the number of basepaired bases
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
    data+=str(pairs['GU']/pair_count)+"\n" #Adds percent of GU pairings

    with open(outFile, 'a') as f:
        f.write(data)