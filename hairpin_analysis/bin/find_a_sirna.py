import argparse
import RNA
import os
import pandas as pd
import re

ops_set = {"A":"U", "U":"A", "C":"G", "G":"C"}
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
        return False
    else:
        return True
    

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file.")
parser.add_argument("-S", "--Sense", help="Whether the siRNA should target the 'Plus' or 'Minus' strand.", default="Plus", choices=["Plus", "Minus"])
parser.add_argument('-l', "--Length_of_oligo", help="Number of bases to parse", default=21)
parser.add_argument("-f", "--foldback_check", help="If flagged, the siRNA is folded using ViennaRNA to see if it forms a hairpin. If it does, the sequence is not used.", default=False, action="store_true")
parser.add_argument("-G", "--exclude_G_stretches", help="If flagged, siRNA with 3 G's in a row are excluded.", default=False, action="store_true")
parser.add_argument("-C", "--exclude_C_stretches", help="If flagged, siRNA with 3 C's in a row are excluded.", default=False, action="store_true")
parser.add_argument("-A", "--exclude_A_stretches", help="If flagged, siRNA with 3 A's in a row are excluded.", default=False, action="store_true")
parser.add_argument("-U", "--exclude_U_stretches", help="If flagged, siRNA with 3 U's in a row are excluded.", default=False, action="store_true")
parser.add_argument("-e", "--ends_with_au", help = "Whether or not the end/start of the siRNA MUST begin with A/U", default=False, action="store_true")
parser.add_argument("-o", "--out", help="Name of output file", default="sirna.csv")
parser.parse_args()
args = parser.parse_args()

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

if runIt:
    if os.path.isfile(outFile):
        data = ""
    else:
        data = "Name,Position,Sense,Sequence,Structure,percentPaired,positionPercent,gcPercent,siRNA\n"
        with open(outFile, 'w') as f:
            f.write(data)
            data = ""

    for i in fastas.keys():
        hit=False

        print("Generating hairpins of sequence: "+i[1:])
        theSeq = fastas[i].replace("T", "U")
        if args.Sense == "Minus":
            theSeq=revc(theSeq)

        print("Folding...")
        fc = RNA.fold_compound(theSeq)
        fc.pf()
        structure=fc.mfe()[0]
        counter=0
        s_count=0
        print("Generating siRNA...")
        while counter+og_length < len(theSeq)+1:
            testSI=theSeq[counter:counter+og_length]

            if args.foldback_check:
                store_it=fold_and_find(testSI)
            else:
                store_it=True

            if args.ends_with_au and ((testSI.endswith("G") or testSI.endswith("C")) or (testSI.startswith("G") or testSI.startswith("C"))) and store_it:
                store_it=False
            
            if args.exclude_G_stretches and store_it and re.search("GGG", testSI):
                store_it=False
                
            if args.exclude_C_stretches and store_it and re.search("CCC", testSI):
                store_it=False

            if args.exclude_A_stretches and store_it and re.search("AAA", testSI):
                store_it=False

            if args.exclude_U_stretches and store_it and re.search("UUU", testSI):
                store_it=False

            if store_it:
                if not hit:
                    hit=True
                s_count+=1
                data=""
                data=i[1:].replace(",", "")+","
                data+=str(counter+1)+","
                data+=args.Sense+","
                data+=testSI.replace("U", "T")+","
                data+=structure[counter:counter+og_length]+","
                data+=str(1-(round(structure[counter:counter+og_length].count(".")/og_length, 3)))+","
                data+=str((counter+1)/len(theSeq))+","
                data+=str(round((testSI.count("G")+testSI.count("C"))/og_length, 3))+","
                data+=revc(theSeq[counter:counter+og_length]).replace("U", "T")
                data+="\n"
                with open(outFile, 'a') as f:
                    f.write(data)
            counter+=1
        if not hit:
            print("No siRNA found with the given parameters :-(")
        else:
            print("Identified "+str(s_count)+" potential "+str(og_length)+"-base siRNA.\n")

    x=pd.read_csv(outFile)
    x=x.sort_values(by=['Name','percentPaired', 'positionPercent', 'gcPercent'])
    x.to_csv(outFile, index=False)
    print("Potential siRNA generated and stored in "+outFile)