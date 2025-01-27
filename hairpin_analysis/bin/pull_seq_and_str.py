import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Path to the input CSV file.") #
parser.parse_args()
args = parser.parse_args()

if args.input:
    try:
        with open(args.input, mode = 'r') as file:
            csvFile = csv.reader(file)
            counter = 0
            for i in csvFile:
                if counter > 0:
                    seq_name = i[0].replace(" ", "-")
                    seq_name = seq_name.replace(",", "-")+"_"+str(counter)
                    new_seq = seq_name+"\n"+i[9]+"\n"+i[10]
                    with open(seq_name.replace(">", ""), 'w') as f:
                        f.write(new_seq)
                counter+=1
    except:
        print("Unable to pull data. Check the input file.")