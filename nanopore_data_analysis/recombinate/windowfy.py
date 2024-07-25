import os
import argparse

#Size of the sliding window to parse each read into
windowLength = 50
#Percent of overlap for the sliding window (0-1). Don't put 1.
overlapPercent = 0

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--reads", help = "Path to the target reads file")
parser.add_argument("-w", "--windowLength", help = "Length of window (default 50)")
parser.add_argument("-o", "--overlap", help = "Percent of overlap for the sliding window")
parser.parse_args()

print("\n")
overlapStart=0

# Read arguments from command line
args = parser.parse_args()
if args.reads:
    print("Parsing reads in % s..." % args.reads)
    try:
        opened_file = open(args.reads, "r")
        reads = opened_file.read()
        opened_file.close()
        reads = reads.strip()
        reads = reads.split(">")
    except:
        print("Unable to open reads file.")
if args.windowLength:
    windowLength = int(args.windowLength)
    print("Splitting reads into windows of % s bases" % windowLength)
if args.overlap:
    overlapStart = int(args.overlap)

try:
    os.mkdir("./recombo")
except:
    print("Recombo direcetory already present...")

if 'reads' in globals():
    print("Beginning reads windowing...\n")
    for read in reads:
        read = read.strip()
        if read != "":
            theRead = read.strip().split("\n")
            readID = theRead[0]
            try:
                os.mkdir("./recombo/"+readID)
            except:
                print("File for read % s already present..." % readID)
            readSEQ = theRead[1]
            readFASTA = ""
            realRange=1+(len(readSEQ)//windowLength)
            if overlapStart > 0:
                extraBases=(realRange-1)*overlapStart
                realRange=1+((len(readSEQ)+extraBases)//windowLength)
            for nth in range(0, realRange):
                if nth < 1:
                    overlap=0
                else:
                    overlap=overlapStart*nth
                windowStart = (nth*windowLength)-overlap
                windowEnd = ((nth+1)*windowLength)-overlap
                if windowEnd > len(readSEQ):
                    windowEnd = len(readSEQ)+1
                readFASTA = readFASTA.strip()+"\n>"+readID+"_"+str(nth)+"\n"+readSEQ[windowStart:windowEnd]

            with open("./recombo/"+readID+'/windowed.fasta', 'w') as f:
                f.write(readFASTA)
else:
    print("Unable to parse reads and/or references")
