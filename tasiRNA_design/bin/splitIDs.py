import argparse
import pandas as pd
# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sirnaFile", help = "Path to sequencing siRNA file", default="sirna.csv")
parser.add_argument("-o", "--outFile", help = "Path to output siRNA fasta file", default="out.fasta")
parser.parse_args()
args = parser.parse_args()
runIt=True

if args.sirnaFile:
    try:
        x=pd.read_csv(args.sirnaFile)
    except:
        print("Unable to load given siRNA file. Abortin.")
        runIt=False

if runIt:
    for i in x.itertuples():
        theList=">"+i[1].replace(" ", "_")+"_"+str(i[2])+"\n"+i[10]+"\n"
        with open(args.outFile, 'a') as f:
            f.write(theList)