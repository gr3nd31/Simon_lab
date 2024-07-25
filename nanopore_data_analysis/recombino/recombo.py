import json
import os
from random import sample
import argparse


# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "Path to the windowed FASTA file")
parser.add_argument("-t", "--tag", help = "Prefix string to the output file")
parser.parse_args()

tag=""

#print("\n")
args = parser.parse_args()
#print(args.input)
if args.tag:
    tag=args.tag
if args.input:
    try:
        opened_file = open(args.input, "r")
        window = opened_file.read()
        opened_file.close()
        new_name = args.input
        new_name = new_name.split(".json")[0]+tag+".csv"
        #print(new_name)
    except:
        print("Input file not found")


final_db="read_id,window,window_start,reference,read_start,read_end,hit_start,hit_end,hit_strand,score,rel_score\n"
windowSet = False

#print("Running recombination check...")
window = json.loads(window)
helios=""
window=window['BlastOutput2']
for i in window:
    alignment = i["report"]["results"]["search"]["hits"]
    read_id = i["report"]["results"]["search"]["query_title"]
    #print(read_id)
    if not windowSet:
        windowLength = i["report"]["results"]["search"]["query_len"]
        windowSet = True
    read_start = (1+int(read_id.split("_")[1]))*windowLength
    if len(alignment) > 0:
        for j in alignment:
            hit_ref = j["description"][0]["title"]
            #print(hit_ref)
            hsps = j["hsps"][0]
            helios=helios+read_id.split("_")[0]+","+str(int(read_id.split("_")[1])+1)+","+str(read_start)+","+hit_ref+","+str(hsps["query_from"])+","+str(hsps["query_to"])+","+str(hsps["hit_from"])+","+str(hsps["hit_to"])+","+str(hsps["hit_strand"])+","+str(hsps["score"])+",0\n"
    else:
        helios=helios+read_id.split("_")[0]+","+str(int(read_id.split("_")[1])+1)+","+str(read_start)+",None,0,0,0,0,0,0,0\n"
final_db+=helios
with open(new_name, 'w') as f:
    f.write(final_db)

