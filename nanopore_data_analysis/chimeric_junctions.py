import pandas
from random import sample
import argparse

virus_file = "Nope"
host_file = "Nope"
reads_file = "Nope"
nn_type = "dna"
snip=10
random_sample=False
sampling_number=2283
save_locale="./"
prefixed=""

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--virus", help = "Path to virus reads_df file.")
parser.add_argument("-c", "--chimeric_analysis", help = "Path to the hits_list file.")
parser.add_argument("-r", "--reads", help = "Path to the fasta file containing the reads")
parser.add_argument("-n", "--nucleic_type", help = "Report as DNA ('dna') or RNA ('rna') sequences.")
parser.add_argument("-l", "--length", help = "Length of sequence to pull on either side of the junction")
parser.add_argument("-o", "--output", help = "Path to output. Default is current directory")
parser.add_argument("-t", "--tag", help = "Prefix tag for files")
parser.parse_args()

print("\n")

# Read arguments from command line
args = parser.parse_args()
if args.virus:
    virus_file = args.virus
if args.chimeric_analysis:
    host_file = args.chimeric_analysis
if args.reads:
    reads_file = args.reads
if args.nucleic_type:
    nn_type = args.nucleic_type
if args.length:
    snip = int(args.length)
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
        #print(target_ids)

        final_file = "read_id,length,viral_start,viral_end,host_start,host_end,orientation,five_prime,junction,three_prime,read\n"

        print("Pulling reads from: "+reads_file)
        reads_dict = {}
        opened_file = open(reads_file, "r")
        fasta_reads = opened_file.read()
        opened_file.close()
        fasta_reads = fasta_reads.split(">")

        for i in fasta_reads:
            if len(i.split(" ")) > 1:
                pick = i.split("\n")[0].split(" ")[0]
            else:
                pick = i.split("\n")[0]
            if pick in target_ids:
                reads_dict[pick] = i.split("\n")[1]

        for i in target_ids:
            #print(i)
            my_read = reads_dict[i]
            my_read = my_read.upper()
            if nn_type == "dna":
                my_read = my_read.replace("U","T")
            else:
                my_read = my_read.replace("T","U")

            virus_align = virus_db[virus_db["read_id"] == i]
            if len(virus_align.index) > 1:
                virus_align = virus_align[virus_align["aligned_length"] == max(virus_align["aligned_length"])]
            if len(virus_align.index) > 1:
                print("Multiple viral alignments detected on read: "+i)
            host_align = host_db[host_db["read_id"] ==i]
            if len(host_align.index) > 1:
                host_align = host_align[host_align["identity"] == max(host_align["identity"])]
                host_align = host_align[host_align["length"] == max(host_align["length"])]
                host_align = host_align[host_align["gaps"] == min(host_align["gaps"])]
                host_align = host_align[host_align["q_start"] == min(host_align["q_start"])]
            if len(host_align.index) > 1:
                host_align = host_align.head(1)
                #print(host_align)
                #print(len(host_align.index))
                print("Multiple host alignments detected on read: "+i)

            #print(len(virus_align.index))
            #print(len(host_align.index))

            final_file = final_file + i + "," + str(len(reads_dict[i])) + "," +str(list(virus_align["q_start"])[0])+","+str(list(virus_align["q_end"])[0])+","+str(list(host_align["q_start"])[0])+","+str(list(host_align["q_end"])[0])+","

            if list(virus_align["q_start"])[0] < list(host_align["q_start"])[0]:
                final_file = final_file+"5prime"
                left_border = list(virus_align["q_end"])[0]
                right_border = list(host_align["q_start"])[0]
            else:
                final_file = final_file+"3prime"
                right_border = list(virus_align["q_start"])[0]
                left_border = list(host_align["q_end"])[0]

            first_prime = my_read[left_border-(snip+1):left_border]
            second_prime = my_read[right_border-1:right_border+snip-1]
            the_junction = my_read[left_border:right_border-1]
            final_file = final_file+","+first_prime+","+the_junction+","+second_prime+","+my_read+"\n"


        with open(save_locale+prefixed+'junction_analysis.csv', 'w') as f:
            f.write(final_file)
