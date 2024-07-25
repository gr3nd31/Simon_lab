import json
import glob
from random import sample
import argparse

#Size of the sliding window to parse each read into
windowLength = 50
#Percent of overlap for the sliding window (0-1). Don't put 1.
overlapPercent = 0

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", help = "Path to the reference file")
parser.add_argument("-r", "--reads", help = "Path to the target reads file")
parser.add_argument("-w", "--windowLength", help = "Length of window (default 50)")
parser.add_argument("-o", "--overlapPercent", help = "Percent of overlap for the sliding window (0-1). Don't put 1.")
parser.add_argument("-t", "--tag", help = "Prefix string to the output file")
parser.parse_args()

print("\n")

# Read arguments from command line
args = parser.parse_args()
if args.genome:
    try:
        opened_file = open(args.genome, "r")
        refs = opened_file.read()
        opened_file.close()
        refs = refs.split(">")
    except:
        print("Unable to open reference file.")
if args.reads:
    try:
        opened_file = open(args.reads, "r")
        reads = opened_file.read()
        opened_file.close()
        reads = reasd.split(">")
    except:
        print("Unable to open reads file.")
if args.windowLength:
    windowLength = args.windowLength
    print("Splitting reads into windows of % s bases" % windowLength)
if args.overlapPercent:
    overlapPercent = args.overlapPercent
    print("Splitting reads into windows of % s bases" % overlapPercent)
if args.tag:
    prefixed=args.tag
    print("Addding prefix of % s to files" % prefixed)

if 'refs' in globals() and 'reads' in globals():
    print("Beginning reads windowing...\n")
    for read in reads:
        theRead = read.split("\n")
        readID = theRead[0]
        readSEQ = theRead[1]
        readFASTA = ""
        for nth in range(1, 1+len(readSEQ)%windowLength):
            print(nth)
else:
    print("Unable to parse reads and/or references")