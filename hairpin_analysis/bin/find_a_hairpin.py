import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sequence", help = "Path to the sequence file. If blank, a random hairpin is generated") #
parser.add_argument("-d", "--dotBracket", help="Name of dotbracketfile") #
parser.add_argument('-m', '--minLength', help="Minimum length a structure should be. Recommended is 11.")
parser.add_argument('-f', '--find', help="Regular expression for the structure to be found. Default is a simple hairpin '\\([\\(\\.]+\\.+[\\)\\.]+\\)'")
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

if args.minLength:
    minLength = int(args.minLength)
else:
    minLength = 0

if args.dotBracket:
    try:
        str_file = open(args.dotBracket, "r")
        structure = str_file.read().strip()
        structure = structure.replace(" ", "")
        str_file.close()
    except:
        print("Unable to read structure file. Aborting.")
        runIt = False
        sequence = ""
else:
    runIt = False
    sequence = ""
    print("Structure file not given. Aborting.")

if args.find:
    try:
        target = args.find
        if target == "":
            print("Target structure not given. Defaulting to a hairpin.")
            target = "\\(+[\\(\\.]+\\.+[\\)\\.]+\\)+"
    except:
        print("Unable to parse target structure. Defaulting to a hairpin.")
        target = "\\(+[\\(\\.]+\\.+[\\)\\.]+\\)+"
else:
    print("Target structure not given. Defaulting to a hairpin.")
    target = "\\(+[\\(\\.]+\\.+[\\)\\.]+\\)+"

if args.prefix:
    prefix = args.prefix
else:
    prefix = "sequence"

if args.out:
    outFile = args.out
else:
    outFile = "sequences.fasta"

if len(sequence) != len(structure):
    print("Sequence and structure are not the saem length. Aborting.")
    runIt = False

if runIt:
    #print(sequence)
    #print(structure)
    finalOut = ""
    p = re.compile(target)
    for m in p.finditer(structure):
        if len(m.group()) >= minLength:
            start = m.start()
            end = m.end()
            finalOut+=">"+prefix+"_"+str(start)+"_"+str(end)+"\n"
            finalOut+=sequence[start:end]+"\n"
    if len(finalOut) > 0:
        with open(outFile, 'a') as f:
            f.write(finalOut)
    else:
        print("Target structure not found. Try again?")