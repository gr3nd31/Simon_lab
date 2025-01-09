import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Path to the input CSV file.") #
parser.add_argument("-s", "--sequence", help = "Name of the sequence file.") #
parser.add_argument("-d", "--dotbracket", help = "Name of the dotbracket file.") #
parser.parse_args()
args = parser.parse_args()

if args.sequence:
    seqFile = args.sequence
else:
    seqFile = "sequence"

if args.dotbracket:
    strFile = args.dotbracket
else:
    strFile = "sequence"

with open(args.input, mode = 'r') as file:
    csvFile = csv.reader(file)
    counter = 0
    for i in csvFile:
        if counter > 0:
            new_seq = seqFile+"_"+str(counter)+".seq"
            new_str = strFile+"_"+str(counter)+".str"
            with open(new_seq, 'w') as f:
                f.write(i[9])
            with open(new_str, 'w') as f:
                f.write(i[10])
        counter+=1