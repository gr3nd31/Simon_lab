# like granulate but only the reads_df function for alignment that is to more than *just* a single fsta (a la CY1)

import json
import glob
from random import sample
import argparse

alignment_file = "Nope"
target_name = "None"
evalue_threshold=1
random_sample=False
sampling_number=2283
save_locale="./"
prefixed=""
max_hsps=10

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--reference", help = "Name of aligned reference to specifically pull reads from.")
parser.add_argument("-a", "--alignment", help = "Path to the alignment json file")
parser.add_argument("-s", "--sampling", help = "Number of alignments to randomly sample")
parser.add_argument("-e", "--evalue_threshold", help = "Minimum evalue score for alignments to be included")
parser.add_argument("-o", "--output", help = "Path to output. Default is current directory")
parser.add_argument("-t", "--tag", help = "Prefix tag for files")
parser.parse_args()

print("\n")

# Read arguments from command line
args = parser.parse_args()
if args.genome:
    genome_file = args.genome
    #print("Path to genome: % s" % args.genome)
if args.alignment:
    alignment_file = args.alignment
    #print("Path to JSON alignment: % s" % args.alignment)
if args.sampling:
    random_sample = True
    sampling_number = int(args.sampling)
    print("Randomly sampling number: % s" % args.sampling)
if args.evalue_threshold:
    evalue_threshold = float(args.evalue_threshold)
    print("Maximum evalue threshold: % s" % args.evalue_threshold)
if args.output:
    if args.output.endswith("/"):
        save_locale = args.output
    else:
        save_locale = args.output+"/"
    print("Saving output files in: % s" % save_locale)
if args.tag:
    prefixed=args.tag
    print("Addding prefix of % s to files" % prefixed)


if alignment_file == "Nope" or "json" not in alignment_file:
    print("JSON alignment file not given.")

else:
    try:
        print("Pulling reads from: "+alignment_file)
        opened_file = open(alignment_file, "r")
        json_file = opened_file.read()
        opened_file.close()
        records= json.loads(json_file)

        aligns_db = "reference,accession,hit_number,read_id,read_length,hit_length,read_hits,hsps_count,hsps_num,q_strand,q_start,q_end,h_strand,h_start,h_end,length,gaps,identity,hseq\n"
        aligned_ids = ""

        records=records['BlastOutput2']
        total_reads=len(records)
        unaligned_number = 0
        aligned_reads = []
        for i in records:
            try:
                if i["report"]["results"]["search"]["message"] != "No hits found":
                    unaligned_number+=1
            except:
                aligned_reads.append(i)

        print("Found "+str(len(aligned_reads))+" reads that aligned.")

        if random_sample:
            print("Randomly sampling "+str(sampling_number)+" out of "+str(aligned_reads)+"...")
            aligned_reads = sample(aligned_reads, sampling_number)
        print("Pulling and sorting IDs from "+str(len(aligned_reads))+" alignments.")
        #print(aligned_reads[0]["report"]["results"]["search"]["hits"][0])
        for i in aligned_reads:
            aligned_number = 1
            aligned_length = 0
            record_it = False
            read_name = i["report"]["results"]["search"]["query_title"]
            read_hits = len(i["report"]["results"]["search"]["hits"])
            for j in i["report"]["results"]["search"]["hits"]:
                #print(read_name+": "+str(len(j["description"])))
                if target_name == "None" or j["description"][0]["title"] == target_name:
                    hsps_count = 1
                    record_it = True
                    for k in j["hsps"]:
                        if k["evalue"] < evalue_threshold:
                            aligns_db = aligns_db+j["description"][0]["title"]+","+j["description"][0]["accession"]+","+str(aligned_number)+","+read_name+","+str(i["report"]["results"]["search"]["query_len"])+","+str(j["len"])+","+str(read_hits)+","+str(len(j["hsps"]))+","+str(hsps_count)+","+k["query_strand"]+","+str(k["query_from"])+","+str(k["query_to"])+","+k["hit_strand"]+","+str(k["hit_from"])+","+str(k["hit_to"])+","+str(k["align_len"])+","+str(k["gaps"])+","+str(k["identity"])+","+k["hseq"]+"\n"
                            hsps_count+=1
                    aligned_number+=1
            if record_it:
                aligned_ids = aligned_ids+i["report"]["results"]["search"]["query_title"]+"\n"

        with open(save_locale+prefixed+'id_list.txt', 'w') as f:
            f.write(aligned_ids)

        with open(save_locale+prefixed+'hits_list.csv', 'w') as f:
            f.write(aligns_db)
    except:
        print("Could not locate file: "+alignment_file)
