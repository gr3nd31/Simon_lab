import argparse
import RNA
import os
import pandas as pd
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

def revc(seq):
    trim=""
    for j in seq:
        trim=ops_set[j]+trim
    return trim

def fold_and_find(seq):
    fc = RNA.fold_compound(seq)
    fc.pf()
    structure=fc.mfe()[0]
    structure = structure.replace(".", "U")
    structure = structure.replace("(", "C")
    structure = structure.replace(")", "G")
    if re.search("C+C+[CU]+C+U+G+[GU]+G+G+", structure):
        #print(structure)
        return False
    else:
        return True
    

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file.") #
parser.add_argument("-S", "--Sense", help="Whether the siRNA should target the 'Plus' or 'Minus' strand.", default="Plus", choices=["Plus", "Minus"])
parser.add_argument('-l', "--Length_of_oligo", help="Number of bases to parse", default=21)
parser.add_argument("-f", "--foldback_check", help="If flagged, the siRNA is folded using ViennaRNA to see if it forms a hairpin. If it does, the sequence is not used.", default=False, action="store_true")
parser.add_argument("-a", "--au_only", help = "Whether or not the end/start of the siRNA MUST begin with A/U", default=False, action="store_true")
parser.add_argument("-t", "--thermoWeight", help="Weighting factor for the percentUnpaired (0-1)", default=0.9)
parser.add_argument("-g", "--gcWeight", help="Weighting factor for the GC content (0-1)", default=0.05)
parser.add_argument("-p", "--positionWeight", help="Weighting factor for the position", default=0.5)
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()
lcounter = 0
runIt=True

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

if args.Length_of_oligo:
    try:
        og_length=int(args.Length_of_oligo)
    except:
        print("Unable to parse given. Is it an integer?")
        runIt=False

# Names the output file
if args.out:
    outFile = args.out
    if not outFile.endswith(".csv"):
        outFile+=".csv"
else:
    outFile = "data.csv"

ops_set = {"A":"U", "U":"A", "C":"G", "G":"C"}

if runIt:
    if os.path.isfile(outFile):
        data = ""
    else:
        data = "Name,Sense,Sequence,Structure,percentPaired,positionPercent,gcPercent,siRNA\n"
        with open(outFile, 'w') as f:
            f.write(data)
            data = ""

    for i in fastas.keys():
        print("Generating hairpins of sequence: "+i[1:])
        theSeq = fastas[i].replace("T", "U")
        if args.Sense == "Minus":
            theSeq=revc(theSeq)

        fc = RNA.fold_compound(theSeq)
        fc.pf()
        structure=fc.mfe()[0]
        counter=0
        while counter+og_length < len(theSeq)+1:
            data=""
            if args.au_only and (theSeq[counter:counter+og_length].endswith("U") or theSeq[counter:counter+og_length].endswith("A")) and (theSeq[counter:counter+og_length].startswith("U") or theSeq[counter:counter+og_length].startswith("A")):
                data=i[1:]+"_"+str(counter+1)+","
                data+=args.Sense+","
                think=theSeq[counter:counter+og_length].replace("U", "T")
                data+=think+","
                data+=structure[counter:counter+og_length]+","
                data+=str(1-(round(structure[counter:counter+og_length].count(".")/og_length, 3)))+","
                data+=str((counter+1)/len(theSeq))+","
                data+=str(round((think.count("G")+think.count("C"))/og_length, 3))+","
                data+=revc(theSeq[counter:counter+og_length]).replace("U", "T")
                data+="\n"
                if args.foldback_check:
                    write_it=fold_and_find(think)
                else:
                    write_it=True
                if write_it:
                    with open(outFile, 'a') as f:
                        f.write(data)
            elif not args.au_only:
                data=i[1:]+"_"+str(counter+1)+","
                data+=args.Sense+","
                think=theSeq[counter:counter+og_length].replace("U", "T")
                data+=think+","
                data+=structure[counter:counter+og_length]+","
                data+=str(1-(round(structure[counter:counter+og_length].count(".")/og_length, 3)))+","
                data+=str((counter+1)/len(theSeq))+","
                data+=str(round((think.count("G")+think.count("C"))/og_length, 3))+","
                data+=revc(theSeq[counter:counter+og_length]).replace("U", "T")
                data+="\n"
                if args.foldback_check:
                    write_it=fold_and_find(think)
                else:
                    write_it=True
                if write_it:
                    with open(outFile, 'a') as f:
                        f.write(data)
            counter+=1

    x=pd.read_csv(outFile)
    x=x.sort_values(by=['percentPaired', 'positionPercent', 'gcPercent'])
    x.to_csv(outFile, index=False)
    print("\nHairpins generated and stored in "+outFile)