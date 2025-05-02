# imports packages
import json
import glob
from random import sample
import argparse

# Set default values of necessary arguments
genome_file = "CY1.fasta"
del_threshold = 7
locale="both.json"
evalue_threshold=1
random_sample=False
sampling_number=2283
save_locale="./"
prefixed=""
max_hsps=10
target_title=""
runit=True

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", help = "Path to the reference file")
parser.add_argument("-a", "--alignment", help = "Path to the alignment json file")
parser.add_argument("-s", "--sampling", help = "Number of alignments to randomly sample")
parser.add_argument("-d", "--deletion_threshold", help = "Minimum deletion size to be recorded. Default is 6")
parser.add_argument("-e", "--evalue_threshold", help = "Minimum evalue score for alignments to be included. Default is 1.")
parser.add_argument("-o", "--output", help = "Path to output. Default is current directory")
parser.add_argument("-t", "--tag", help = "Prefix tag for files")
parser.add_argument("-r", "--reference", help = "Name to put in the 'reference' column of the reads_df.csv. Default is name of reference as it is presented in the json description.")
parser.parse_args()
print("\n")

# Read arguments from command line
args = parser.parse_args()
# Gets the genome fasta file
if args.genome:
    genome_file = args.genome
# Gets the aligment json file
if args.alignment:
    locale = args.alignment
#Sets random sampling to TRUE and pulls the number
if args.sampling:
    try:
        random_sample = True
        sampling_number = int(args.sampling)
        print("Randomly sampling number: % s" % args.sampling)
    except:
        random_sample = False
        print("Unable to read integer value for random sampling. Defaulting to the full dataset")
# Sets the minimum number of bases to call something a deleted sequence
if args.deletion_threshold:
    del_threshold = int(args.deletion_threshold)
    print("Minimum deletion length: % s" % args.deletion_threshold)
# Sets a default evalue threshold
if args.evalue_threshold:
    try:
        evalue_threshold = float(args.evalue_threshold)
        print("Maximum evalue threshold: % s" % args.evalue_threshold)
    except:
        print("Unable the read e-value threshold")
# Sets output file location
if args.output:
    if args.output.endswith("/"):
        save_locale = args.output
    else:
        save_locale = args.output+"/"
    print("Saving output files in: % s" % save_locale)
# Adds a prefix to the file
if args.tag:
    prefixed=args.tag
    print("Addding prefix of % s to files" % prefixed)
# Sets a reference name to specifically pull from the alignment
if args.reference:
    target_title=args.reference
else:
    target_title="none"


print("Comparing reads in "+locale+" with "+genome_file+"...\n")

try:
    #Opens the genome file
    trick = open(genome_file, "r")
    genome = trick.read()
    genome2=genome.strip().split("\n")
    trick.close()
    # Creates second genome copy for detecting insertions
    trickk=""
    for i in range(1,len(genome2)):
        trickk+=genome2[i]
    genome=trickk

    # Creates a genome dictionary to track various errors
    genome_dict = {}
    start=1
    for i in genome:
#                               Letter, hit, hit_letter,                            mismatch, deletion, insert, insert_letter
        genome_dict[str(start)] = [i.upper(), 0, {"G":0, "A":0, "T":0, "C":0, "U":0, "E":0}, 0, 0, 0, {"G":0, "A":0, "T":0, "C":0, "U":0, "E":0}]
        start += 1
except:
    print("Unable to open the genome file. Aborting.")
    runit=False

#Creates headers for the output csv files
del_df="Number,read_id,orientation,sense,position,size,del_seq\n"
ins_df="Number,read_id,orientation,sense,position,size,ins_seq\n"
read_df="Number,read_id,reference,hsps,bit_score,evalue,coverage,q_strand,h_strand,orientation,read_length,length,h_start,h_end,q_start,q_end,skipped,deletion_num,mismatch_num,insert_num\n"

if runit:
    print("Parsing file: "+locale)
    # opens the alignment json file
    try:
        opened_file = open(locale, "r")
        json_file = opened_file.read()
        opened_file.close()
        records= json.loads(json_file)
    except:
        print("Unable to open the alignment file. Aborting.")
        runit=False

    # Initiates the aligned read, unaligned reads, and error counts
    unaligned_reads=0
    aligned_reads=0
    errors=0

if runit:
    try:
        # Writes the reads and deletions files, since these are appended. Since the genome file stays small, it does not need to be appended.
        with open(save_locale+prefixed+'reads_df.csv', 'w') as f:
            f.write(read_df)

        with open(save_locale+prefixed+'del_df.csv', 'w') as f:
            f.write(del_df)
    except:
        print("Unable to save files at given location. Aborting.")
        runit=False

if runit:
    try:
        # Accesses total alignments
        records=records['BlastOutput2']
        total_reads=len(records)
        # Randomly samples reads, if random sampling is on.
        # Note: If a large amount of reads DO NOT align to the reference genome, sampling may result in fewer reads than desired.
        if random_sample:
            records=sample(records, sampling_number)
            print("Randomly sampling "+str(sampling_number)+" out of "+str(total_reads)+"...")
        print("Found a total of "+str(total_reads)+" reads in blast file")

        # Initiates a tracker to alert the user to the current analysis count
        spot_read=0
        #Iterates through the sample reads
        for i in records:
            j=i["report"]["results"]["search"]
            try:
                # Adds to the read counter
                spot_read+=1
                # Alerts the user ever 10,000 reads
                if spot_read%10000 == 0:
                    print("Analyzing read: "+str(spot_read))
                #initiates an hsps number to record the number of hsps's per read
                hsps_number=0
                #Records the read name
                read_id = j["query_title"]
                # Records the length of the read
                read_length = j["query_len"]
                if target_title=="none":
                    target_title=j["hits"][0]["description"][0]["title"]
                # iterates through each hsps alignment for the read
                read_data_raw = j["hits"][0]["hsps"]
                for hsps in read_data_raw:
                    read_data=hsps
                    # counts the number of HSPS
                    hsps_number+=1
                    # intiates a 'skip' variable to skip reads that do not align to anything or have an evalue score above the given threshold
                    skipit="False"
                    try:
                        if read_data['evalue'] < evalue_threshold:
                            skipit="False"
                    except:
                        if "e-" in read_data["evalue"]:
                            skipit="False"
                        
                    #records the starting hit position
                    h_from=read_data['hit_from']
                    #records the ending hit position
                    h_to=read_data['hit_to']
                    #records the starting read position
                    q_from=read_data['query_from']
                    #records the ending read position
                    q_to=read_data['query_to']
                    #records the orientation of the read
                    q_strand=read_data['query_strand']
                    #records the orientation of the aligment
                    h_strand=read_data['hit_strand']
                    #records the alignment length
                    span=read_data['align_len']

                    # Initates the inserts, mismatches, and deletions for the hsps
                    current_inserts = 0
                    current_misses = 0
                    current_dels = 0

                    #Picks the longest alignment if two hsps's overlap
                    if len(read_data_raw) > 1 and hsps_number > 1:
                        old_h_from = read_data_raw[hsps_number-2]['hit_from']
                        old_h_to = read_data_raw[hsps_number-2]['hit_to']

                    if h_from > h_to:
                        direction="Reverse"
                        ticker = -1
                    else:
                        direction="Forward"
                        ticker = 1

                    if len(read_data_raw) > 1 and hsps_number > 1:
                        if h_from >= old_h_from and h_from <= old_h_to:
                            skipit="True"
                        elif old_h_to >= h_from and old_h_to <= h_to:
                            skipit="True"

                    # If the read is NOT skipped, the alignment is parsed
                    if skipit == "False":
                        # Pulls the hit sequence and the query sequence
                        hseq=read_data['hseq']
                        qseq=read_data['qseq']

                        # initiates states to detect either a deleted or inserted sequence
                        del_count=1
                        del_state=False
                        insert_detected = False
                        del_seq=""
                        #Makes sure they're the same size
                        if len(hseq) == len(qseq):
                            new_qseq=""
                            #Removes the inserts from the query by seeing if there is an '-' in the hit sequence
                            for k in range(0,len(qseq)):
                                if hseq[k]=="-":
                                    new_qseq+=""
                                else:
                                    new_qseq+=qseq[k]
                            # hseq = hit sequence with no inserts
                            # gymn_seq = hit sequence WITH inserts
                            # new_seq = query sequence with deletions but without inserts
                            # qseq = query sequence with deletions and inserts
                            gymn_seq = hseq
                            hseq=hseq.replace("-", "")
                        try:
                            # Iterates through the query sequence
                            for k in range(0,len(new_qseq)):
                                #Tracks frequency of position occuring
                                the_step = "position"
                                genome_dict[str((ticker*k)+h_from)][1]+=1

                                #Tracks mismatch freqncy and alternate base by finding places the new_seq mismatches with the hseq but isn't a deletion
                                the_step = "mismatch"
                                if hseq[k].upper() != new_qseq[k].upper() and new_qseq[k] != "-" and not del_state:
                                    try:
                                        if genome_dict[str((ticker*k)+h_from)][0] != new_qseq[k].upper():
                                            genome_dict[str((ticker*k)+h_from)][2][new_qseq[k].upper()]+=1
                                            genome_dict[str((ticker*k)+h_from)][3]+=1
                                            current_misses+=1
                                    except:
                                        pass
                                        #print(new_qseq[k].upper())
                                        #print(genome_dict[str((ticker*k)+h_from)][2])
                                    

                                #Tracks insertions and inserted bases
                                the_step = "insert"
                                if gymn_seq[0] == "-" and qseq[k] != "-":
                                    genome_dict[str((ticker*k)+h_from)][5]+=1
                                    try:
                                        genome_dict[str((ticker*k)+h_from)][6][qseq[k].upper()]+=1
                                    except:
                                        pass
                                        #print(qseq[k].upper())
                                    insert_detected=True
                                    current_inserts+=1
                                    if ticker==1:
                                        gymn_seq=gymn_seq[1:len(gymn_seq)]
                                        #qseq=qseq[1:len(qseq)]
                                    else:
                                        gymn_seq=gymn_seq[0:len(gymn_seq)-1]
                                        #qseq=qseq[0:len(qseq)-1]
                                if ticker==1:
                                    gymn_seq=gymn_seq[1:len(gymn_seq)]
                                    #qseq=qseq[1:len(qseq)]
                                else:
                                    gymn_seq=gymn_seq[0:len(gymn_seq)-1]
                                    #qseq=qseq[0:len(qseq)-1]

                                #Tracks deletions
                                # Detects the start of a new deletion
                                the_step = "start_deletions"
                                if new_qseq[k]=="-" and del_state==False:
                                    genome_dict[str((ticker*k)+h_from)][4]+=1
                                    del_state=True
                                    current_dels+=1
                                    del_start=(ticker*k)+h_from
                                    del_seq+=genome_dict[str((ticker*k)+h_from)][0]
                                #Detects a continued deletion
                                elif new_qseq[k]=="-" and del_state==True:
                                    current_dels+=1
                                    genome_dict[str((ticker*k)+h_from)][4]+=1
                                    del_seq+=genome_dict[str((ticker*k)+h_from)][0]
                                # Detects the end of a string of deletions
                                elif new_qseq[k]!="-" and del_state==True:
                                    del_state=False
                                    if hseq[k].upper() != new_qseq[k].upper():
                                        try:
                                            if genome_dict[str((ticker*k)+h_from)][0] != new_qseq[k].upper():
                                                genome_dict[str((ticker*k)+h_from)][2][new_qseq[k].upper()]+=1
                                                genome_dict[str((ticker*k)+h_from)][3]+=1
                                                current_misses+=1
                                        except:
                                            #print(new_qseq[k].upper())
                                            #print(genome_dict[str((ticker*k)+h_from)][2])
                                            pass

                                    if len(del_seq) >= del_threshold:
                                        if (ticker*k)+h_from > len(genome)/2:
                                            del_orient = "Back"
                                        else:
                                            del_orient = "Front"
                                        # Appends the deletions to the existing csv file
                                        with open(save_locale+prefixed+'del_df.csv', 'a') as f:
                                            f.write(str(del_count)+","+read_id+","+direction+","+del_orient+","+str(del_start)+","+str(len(del_seq))+","+del_seq+"\n")
                                        del_seq=""
                                        del_count+=1
                                if insert_detected and gymn_seq[0] != "-":
                                    genome_dict[str((ticker*k)+h_from)][6]["E"]+=1
                                    insert_detected = False
                        except:
                            #print(the_step)
                            errors+=1
                            #print(new_seq)
                        aligned_reads+=1
                    # Appends the read to the existing csv file
                    with open(save_locale+prefixed+'reads_df.csv', 'a') as f:
                        f.write(str(spot_read)+","+read_id+","+target_title+","+str(hsps_number)+","+str(read_data["bit_score"])+","+str(read_data["evalue"])+","+str(100*(span/len(genome)))+","+q_strand+","+h_strand+","+direction+","+str(read_length)+","+str(read_data['align_len'])+","+str(h_from)+","+str(h_to)+","+str(q_from)+","+str(q_to)+","+skipit+","+str(current_dels)+","+str(current_misses)+","+str(current_inserts)+"\n")
                    max_insertion_length = 0
                if random_sample and aligned_reads >= sampling_number:
                    break

            #If no alignment is found, the script passes
            except:
                unaligned_reads+=1
                pass
    except:
        print("Yeah, something went wrong")

    # prints the number and percent of aligned reads
    print("Aligned reads: "+str(aligned_reads)+" ("+str(100*round(aligned_reads/total_reads,2))+"%)")
    # If aligned reads exist, the genome file is written
    if aligned_reads > 0:
        print("Aligned read errors: "+str(errors)+" ("+str(100*round(errors/aligned_reads,2))+"%)")
        print("\n")
        print("Generating genome df...")
        genome_df="Position_num,Position_nt,number_hit,Other_nt,number_mismatched,number_deleted,number_inserted,inserted_nt\n"
        for i in genome_dict:
            genome_df=genome_df+i+","+genome_dict[i][0]+","+str(genome_dict[i][1])+",G:"+str(genome_dict[i][2]["G"])+"_C:"+str(genome_dict[i][2]["C"])+"_A:"+str(genome_dict[i][2]["A"])+"_T:"+str(genome_dict[i][2]["T"])+"_U:"+str(genome_dict[i][2]["U"])+","+str(genome_dict[i][3])+","+str(genome_dict[i][4])+","+str(genome_dict[i][5])+",G:"+str(genome_dict[i][6]["G"])+"_C:"+str(genome_dict[i][6]["C"])+"_A:"+str(genome_dict[i][6]["A"])+"_T:"+str(genome_dict[i][6]["T"])+"_U:"+str(genome_dict[i][6]["U"])+"_E:"+str(genome_dict[i][6]["E"])+"\n"

        print("Writing files...")
        with open(save_locale+prefixed+'genome_df.csv', 'w') as f:
            f.write(genome_df)

    else:
        print("No reads found :-(")
        print("\n")

    print("Granulation complete.\n")
