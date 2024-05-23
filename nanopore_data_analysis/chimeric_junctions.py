import pandas
from random import sample
import argparse

virus_file = "Nope"
host_file = "Nope"
reads_file = "Nope"
random_sample=False
sampling_number=2283
save_locale="./"
prefixed=""

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--virus", help = "Path to virus reads_df file.")
parser.add_argument("-h", "--host", help = "Path to the hits_list file.")
parser.add_argument("-r", "--reads", help = "Path to the fasta file containing the reads")
parser.add_argument("-s", "--sampling", help = "Number of alignments to randomly sample")
parser.add_argument("-o", "--output", help = "Path to output. Default is current directory")
parser.add_argument("-t", "--tag", help = "Prefix tag for files")
parser.parse_args()

print("\n")

# Read arguments from command line
args = parser.parse_args()
if args.virus:
    virus_file = args.virus
if args.host:
    host_file = args.host
if args.sampling:
    random_sample = True
    sampling_number = int(args.sampling)
    print("Randomly sampling number: % s" % args.sampling)
if args.reads:
    reads_file = args.reads
if args.output:
    if args.output.endswith("/"):
        save_locale = args.output
    else:
        save_locale = args.output+"/"
    print("Saving output files in: % s" % save_locale)
if args.tag:
    prefixed=args.tag
    print("Addding prefix of % s to files" % prefixed)


if virus_file == "Nope" or "csv" not in virus_file:
    print("Viral reads_df file not given.")
elif host_file == "Nope" or "csv" not in host_file:
    print("Chimeric hits_list file not given.")
elif reads_file == "Nope" or "fasta" not in reads_file:
    print("Fasta file not goven.")
else:
    #try:
        print("Using viral reads alignment from: "+virus_file)
        virus_db = pandas.read_csv(virus_file)
        virus_read_ids = virus_db["read_id"].unique()

        print("Using chimeric reads alignment from: "+host_file)
        host_db = pandas.read_csv(host_file)
        host_reads_ids = host_db["read_id"].unique()

        target_ids = list(set(virus_read_ids).intersection(host_reads_ids))

        print("Pulling reads from: "+reads_file)
        reads_dict = {}
        opened_file = open(virus_file, "r")
        fasta_reads = opened_file.read()
        opened_file.close()
        fasta_reads = fasta_reads.split(">")
        for i in fasta_reads:
            if i.split("\n")[0] in target_ids:
                reads_dict[i.split("\n")[0]] = i.split("\n")[0]
                print(i.split("\n")[0])


#        with open(save_locale+prefixed+'id_list.txt', 'w') as f:
#            f.write(aligned_ids)

#        with open(save_locale+prefixed+'hits_list.csv', 'w') as f:
#            f.write(aligns_db)
