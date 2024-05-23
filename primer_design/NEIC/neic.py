import random
import argparse

# Default variables
fasta_file = "Nope"
number_of_oligo = 0
length_of_oligo = 0
attempts = 0
out_prefix="neics.txt"

# Real variables
base_list = ["G","A","C","U"]
limited = False

# Define a function
def oligo_invert(oligo):
    zim = ""
    for i in oligo[::-1]:
        if i == "G":
            zim += "C"
        elif i == "C":
            zim += "G"
        elif i == "U":
            zim += "A"
        elif i == "A":
            i += "U"
    return zim


# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help = "Path to fasta file.")
parser.add_argument("-n", "--number_of_oligos", help = "Number of oligos to generate (min of 1).")
parser.add_argument("-l", "--length_of_oligos", help = "Length of bases to generate (min of 2.")
parser.add_argument("-a", "--attempts", help = "Number of attempts before stopping (min of 1).")
parser.add_argument("-o", "--output", help = "Name of output file. Defaults is 'neics.csv'")
parser.parse_args()

# Read arguments from command line
args = parser.parse_args()
if args.fasta:
    fasta_file = args.fasta
if args.number_of_oligos:
    number_of_oligo = int(args.number_of_oligos)
if args.length_of_oligos:
    length_of_oligo = int(args.length_of_oligos)
if args.attempts:
    limited = True
    attempts = args.attempts
if args.output:
    out_prefix = args.output

# Checks that parameters work
if fasta_file == "Nope" or "fasta" not in fasta_file:
    print("Target fasta file not found.")
elif number_of_oligo < 1:
    print("Number of oligos not designated")
elif length_of_oligo < 2:
    print("Length of extensions not given")
else:
    try:
        oligo_list = []
        check_list = []
        if not limited:
            attempts = 4**length_of_oligo

        print("Finding "+str(number_of_oligo)+" oligos with a length of "+str(length_of_oligo)+" do not bind to: "+fasta_file)
        #print("\n")
        opened_file = open(fasta_file, "r")
        fasta_reads = opened_file.read()
        opened_file.close()

        if len(fasta_reads.split("\n")) > 2:
            fasta_seq = ""
            for i in fasta_reads.split("\n")[1:]:
                fasta_seq+=i
        elif len(fasta_reads.split("\n")) == 2:
            fasta_seq = fasta_reads.split("\n")[1]
        #print(fasta_seq)
        fasta_seq = fasta_seq.upper()
        if "T" in fasta_seq:
            print("DNA sequence detected, converting to RNA...")
            fasta_seq = fasta_seq.replace("T", "U")
        number_of_attempts = 0
        while len(oligo_list) < number_of_oligo and number_of_attempts < attempts:
            o_check = False
            possible_oligo = ""
            for i in range(0,length_of_oligo):
                possible_oligo += random.choice(base_list)

            #print(possible_oligo)
            if oligo_invert(possible_oligo) not in fasta_seq and possible_oligo not in oligo_list:
                o_check = True
                mis_list = []
            if o_check:
                oligo_list.append(possible_oligo)
            if possible_oligo not in check_list:
                number_of_attempts += 1
                check_list.append(possible_oligo)

        if number_of_attempts == attempts:
            print("Reached maximum number of attempts ("+str(attempts)+"). Stopping oligo generation")
        #print(len(oligo_list))
        #print(number_of_attempts)
        #print(oligo_list)

        if ".txt" not in out_prefix:
            out_prefix += ".txt"

        trip = ""
        for i in oligo_list:
            trip += i+"\n"

        with open(out_prefix, 'w') as f:
            f.write(trip)
        print("Oligo report completed.\nFound "+str(len(oligo_list))+" oligos with complement sequences that are not found in "+fasta_file+".")
    except:
        print("Something went wrong...")
