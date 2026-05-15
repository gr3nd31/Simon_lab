import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--peCSVFile", help = "Path to the CSV output file containing the PE values.") #
parser.add_argument("-o", "--output", help = "Name of the output directory.", default="PE_files") #
parser.parse_args()
args = parser.parse_args()
runIt=True

try:
    structs=pd.read_csv(args.peCSVFile)
except:
    print("Unable to open the given CSV file. Aborting")
    runIt=False

try:
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    else:
        print("PE output directory already exsists.")
except:
    print("Unable to make output directory. Aborting")

if runIt:
    for i in range(0,len(structs['Name'])):
        fileContents=">"+structs['Name'][i].replace(">", "")+"\n"+structs["Sequence"][i]+"\n"+structs["PE"][i]
        outFile=args.output+"/"+structs['Name'][i].replace(">", "")+".txt"
        with open(outFile, 'a') as f:
            f.write(fileContents)