import argparse
import pandas as pd

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Path to the RNACanvas file")
parser.add_argument("-d", "--data", help = "Path to the CSV file containing data")
parser.add_argument("-p", "--palette", help = "Viridis palette to use. Default is 'inferno'")
parser.add_argument("-o", "--output", help = "Path to output. Default is current directory")
parser.add_argument("-t", "--threshold", help = "Threshold to begin coloring at.")
parser.parse_args()

# Read arguments from command line
args = parser.parse_args()

runIt = True

if args.input:
    try:
        trick = open(args.input, "r")
        canvas = trick.read()
        trick.close()
    except:
        print("Unable to read RNA canvas file. Check the path and try again")
        runIt = False
else:
    print("RNA canvas file not given.")
    runIt = False

if args.data:
    try:
        data=pd.read_csv(args.data)
        data.columns = ['number', 'data', 'hex']
    except:
        print("Unable to read data file. Check the path and try again")
        runIt = False
else:
    print("Data file not given.")
    runIt = False

if args.output:
    outputFile=args.output
    if not outputFile.endswith('rnacanvas'):
        print("Adding the '.rnacanvas' extension to "+args.output)
        outputFile+=".rnacanvas"
elif args.input:
    outputFile="brained_"+args.input

if runIt:
    the_start = ""
    counter = 1
    the_end = canvas.split("</text><line")[1]
    canvas=canvas.split("</text><line")[0]
    if len(canvas.split("</tspan>")) < len(data):
        print("Given data is longer than bases in canvas file. Check files and try again.")
    else:
        for i in canvas.split("</tspan>"):
            if counter in data['number'] and len(i) > 0:
                the_start+=i.split('fill=\\"')[0]
                the_start+='fill=\\"'
                the_start+=data[data['number'] == counter]['hex'].values[0]+''
                the_start+=i.split('fill=\\"')[1][7:]
                the_start+="</tspan>"
            else:
                the_start+=i+"</tspan>"
            counter+=1
            #print(the_start)

    the_start+="</text><line"+the_end            
    with open(outputFile, "w") as f:
            f.write(the_start)