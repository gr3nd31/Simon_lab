import re
import argparse
import RNA

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence and structure file") #
parser.add_argument('-S', "--Structure", help = "Add stucture to the output file.", action = 'store_true')
parser.add_argument('-m', '--minLength', help="Minimum length a structure should be. Recommended is 11.")
parser.add_argument('-b', '--bonds', help="Minimum number of pairs required. Default minimum is 3.")
parser.add_argument('-F', '--Force', help='If flagged, the subsequence is again folded to verify the sequence can fold on its own and includes no additional sequences.',
                     action='store_true')
parser.add_argument('-f', '--find', help="Regular expression for the structure to be found. Default is a simple hairpin 'C+C+[CU]+C+U+G+[GU]+G+G+'. A YSS would be 'C+[CU]+U+G+[GU]+C+[CU]+G+[GU]+G+'")
parser.add_argument('-p', '--prefix', help="Prefix to be added to the sequence.")
parser.add_argument("-o", "--out", help="Name of output file") #
parser.parse_args()
args = parser.parse_args()

runIt = True
# Pulls a target sequence
if args.sequence:
    try:
        # Pulls data from file and splits it
        seq_file = open(args.sequence, "r")
        datum = seq_file.read().strip().split("\n")
        seq_file.close()
        # Pulls sequence data and replaces T's to U's
        sequence = datum[1].upper().replace(" ", "")
        sequence = sequence.replace("T", "U")

        #Pulls structure from file and converts to CUG nomenclature
        structure = datum[2].replace(" ", "")
        structure = structure.replace(".", "U")
        structure = structure.replace("(", "C")
        structure = structure.replace(")", "G")

         # Defines the prefix of the subsequence
        if args.prefix:
            prefix = args.prefix
        else:
            prefix = datum[0].replace(">", "")
            prefix = prefix.replace(",", " -")
            prefix = prefix.replace("_", "-")
        
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
    minLength = 8

# Defines the minimum number of basepairs
if args.bonds:
    try:
        minBonds = int(args.bonds)
    except:
        minBonds = 3
        print("Noninteger given for minimum bonds. Setting to default of 3.")
else:
    minBonds = 3

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

# Defines the name of the out file
if args.out:
    outFile = args.out
else:
    outFile = "sequences.fasta"

if args.Force:
    print("Checking substructues to verify correct extraction.")

#Checks is the len of the structure file is the same as a sequence
if len(sequence) != len(structure):
    print("Sequence and structure are not the same length. Aborting.")
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
            dStart = start
            end = m.end()
            dEnd = end
            # Splits the substructure into the 5prime and 3prime segments
            if args.Force:
                fc = RNA.fold_compound(sequence[start:end])
                fc.pf()
                subStructure = fc.mfe()[0].replace(",", ".")
                subStructure = subStructure.replace(".", "U")
                subStructure = subStructure.replace("(", "C")
                subStructure = subStructure.replace(")", "G")
                for n in p.finditer(subStructure):
                    if len(n.group()) >= minLength and subStructure[n.start():n.end()].count("C") >= minBonds:
                        start = n.start()+dStart
                        end = n.end()+dStart
                        finalOut+=">"+prefix+"_"+str(start+1)+"_"+str(end)+"\n"
                        finalOut+=sequence[start:end]+"\n"
                        if args.Structure:
                            finalOut+=subStructure[n.start():n.end()].replace("U",".").replace("C","(").replace("G",")")+"\n"
                    else:
                        continue

            #If the trimmed sequence is long enough and contains enough basepairs, the substructure sequence is pulled.
            elif len(sequence[start:end]) >= minLength and structure[start:end].count("C") >= minBonds:
                finalOut+=">"+prefix+"_"+str(start+1)+"_"+str(end)+"\n"
                finalOut+=sequence[start:end]+"\n"
                if args.Structure:
                    finalOut+=structure[start:end].replace("U",".").replace("C","(").replace("G",")")+"\n"
    # If structures were found, they are saved
    if len(finalOut) > 0:
        with open(outFile, 'a') as f:
            f.write(finalOut)
    else:
        print("Target structure not found. Try again?")
