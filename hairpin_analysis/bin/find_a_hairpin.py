import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file. If blank, a random hairpin is generated") #
parser.add_argument("-d", "--dotBracket", help="Name of dotbracketfile") #
parser.add_argument('-m', '--minLength', help="Minimum length a structure should be. Recommended is 11.")
parser.add_argument('-b', '--bonds', help="Minimum number of pairs required. Default minimum is 3.")
parser.add_argument('-f', '--find', help="Regular expression for the structure to be found. Default is a simple hairpin 'C+C+[CU]+C+U+G+[GU]+G+G+'")
parser.add_argument('-p', '--prefix', help="Prefix to be added to the sequence.")
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()

runIt = True
# Pulls a target sequence
if args.sequence:
    try:
        seq_file = open(args.sequence, "r")
        sequence = seq_file.read().strip()
        sequence = sequence.upper().replace(" ", "")
        sequence = sequence.replace("T", "U")
        seq_file.close()
    except:
        print("Unable to read sequence file. Aborting.")
        runIt = False
        sequence = ""
else:
    runIt = False
    sequence = ""
    print("Sequence file not given. Aborting.")

# Defines a mininum length of sequence
if args.minLength:
    minLength = int(args.minLength)
else:
    minLength = 0

# Defines the minimum number of basepairs
if args.bonds:
    try:
        minBonds = int(args.bonds)
    except:
        minBonds = 3
        print("Noninteger given for minimum bonds. Setting to default of 3.")
else:
    minBonds = 3

# Loads in the dot bracket structure and replaceswith more regex friendly notation
if args.dotBracket:
    try:
        str_file = open(args.dotBracket, "r")
        structure = str_file.read().strip()
        structure = structure.replace(" ", "")
        structure = structure.replace(".", "U")
        structure = structure.replace("(", "C")
        structure = structure.replace(")", "G")
        str_file.close()
    except:
        print("Unable to read structure file. Aborting.")
        runIt = False
        sequence = ""
else:
    runIt = False
    sequence = ""
    print("Structure file not given. Aborting.")

# Defines the target structure. '.' : 'U', '(' : 'C', ')': 'G'
if args.find:
    try:
        target = args.find
        if target == "":
            print("Target structure not given. Defaulting to a hairpin.")
            target = "C+C+[CU]+C+U+G+[GU]+G+G+"
    except:
        print("Unable to parse target structure. Defaulting to a hairpin.")
        target = "C+C+[CU]+C+U+G+[GU]+G+G+"
else:
    print("Target structure not given. Defaulting to a hairpin.")
    target = "C+C+[CU]+C+U+G+[GU]+G+G+"

# Defines the prefix of the subsequence
if args.prefix:
    prefix = args.prefix
else:
    prefix = "sequence"
# Defines the name of the out file
if args.out:
    outFile = args.out
else:
    outFile = "sequences.fasta"

#Checks is the len of the structure file is the same as a sequence
if len(sequence) != len(structure):
    print("Sequence and structure are not the saem length. Aborting.")
    runIt = False

# Only runs if the necessary parameters are defined
if runIt:
    # Starts the final string
    finalOut = ""
    # Compiles the regex
    p = re.compile(target)
    # Iterates through the found substructures that match the regex
    for m in p.finditer(structure):
        # Checks that the substructure is long enough
        if len(m.group()) >= minLength:
            # Defines the start and end positions to be extracted
            start = m.start()
            end = m.end()
            # Splits the substructure into the 5prime and 3prime segments
            trick = re.split("CU+G", structure[start:end])
            # Proceeds if there is only one apical loop 
            if len(trick) == 2:
                # Pulls the 5prime and 3prime ends and trims any excess
                p5 = trick[0]
                p3 = trick[1]
                if p5.count("C") > p3.count("G"):
                    clock = p3.count("G")+1
                    for i in p5[::-1]:
                        if i == "C":
                            clock-=1
                            if clock <= 0:
                                start+=1
                        elif i == "U" and clock <= 1:
                            start+=1
                elif p5.count("C") < p3.count("G"):
                    clock = p5.count("C")+1
                    for i in p3:
                        if i == "G":
                            clock-=1
                            if clock <= 0:
                                end-=1
                        elif i == "U" and clock <= 1:
                            end-=1
                else:
                    start = m.start()
                    end = m.end()
            #If the trimmed sequence is long enough and contains enough basepairs, the substructure sequence is pulled.
            if len(sequence[start:end]) >= minLength and structure[start:end].count("C") >= minBonds:
                finalOut+=">"+prefix+"_"+str(start)+"_"+str(end)+"\n"
                finalOut+=sequence[start:end]+"\n"
                # Uncomment to also record the structure
                #finalOut+=structure[start:end].replace("U",".").replace("C","(").replace("G",")")+"\n"
    # If structures were found, they are saved
    if len(finalOut) > 0:
        with open(outFile, 'a') as f:
            f.write(finalOut)
    else:
        print("Target structure not found. Try again?")