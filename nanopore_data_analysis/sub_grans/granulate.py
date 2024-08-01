import json
import glob
from random import sample
import argparse

genome_file = "CY1.fasta"
del_threshold = 6
locale="both.json"
evalue_threshold=1
random_sample=False
sampling_number=2283
save_locale="./"
prefixed=""
max_hsps=10
target_title=""

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", help = "Path to the reference file")
parser.add_argument("-a", "--alignment", help = "Path to the alignment json file")
parser.add_argument("-s", "--sampling", help = "Number of alignments to randomly sample")
parser.add_argument("-d", "--deletion_threshold", help = "Minimum deletion size to be recorded")
parser.add_argument("-e", "--evalue_threshold", help = "Minimum evalue score for alignments to be included")
parser.add_argument("-o", "--output", help = "Path to output. Default is current directory")
parser.add_argument("-t", "--tag", help = "Prefix tag for files")
parser.add_argument("-r", "--reference", help = "Name of reference as it is presented in the json description.")
parser.parse_args()

print("\n")

# Read arguments from command line
args = parser.parse_args()
if args.genome:
    genome_file = args.genome
if args.alignment:
    locale = args.alignment
if args.sampling:
    random_sample = True
    sampling_number = int(args.sampling)
    print("Randomly sampling number: % s" % args.sampling)
if args.deletion_threshold:
    del_threshold = int(args.deletion_threshold)
    print("Minimum deletion length: % s" % args.deletion_threshold)
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
if args.reference:
    target_title=args.reference


print("Comparing reads in "+locale+" with "+genome_file+"...\n")

trick = open(genome_file, "r")
genome = trick.read()
genome2=genome.strip().split("\n")
trick.close()

trickk=""
for i in range(1,len(genome2)):
    trickk+=genome2[i]
genome=trickk

called=True

genome_dict = {}
start=1
for i in genome:
#                    Letter, hit, hit_letter, mismatch, deletion, insert, insert_letter
    genome_dict[str(start)] = [i.upper(), 0, {"G":0, "A":0, "T":0, "C":0, "U":0, "E":0}, 0, 0, 0, {"G":0, "A":0, "T":0, "C":0, "U":0, "E":0}]
    start += 1

del_df="Number,read_id,orientation,sense,position,size,del_seq\n"
ins_df="Number,read_id,orientation,sense,position,size,ins_seq\n"
read_df="Number,read_id,hsps,bit_score,evalue,coverage,q_strand,h_strand,orientation,read_length,length,h_start,h_end,q_start,q_end,skipped,deletion_num,mismatch_num,insert_num\n"

with open(save_locale+prefixed+'reads_df.csv', 'w') as f:
    f.write(read_df)

with open(save_locale+prefixed+'del_df.csv', 'w') as f:
    f.write(del_df)

print("Parsing file: "+locale)

opened_file = open(locale, "r")
json_file = opened_file.read()
opened_file.close()
records= json.loads(json_file)

unaligned_reads=0
aligned_reads=0
errors=0

try:
    records=records['BlastOutput2']
    total_reads=len(records)
    if random_sample:
        records=sample(records, sampling_number)
        print("Randomly sampling "+str(sampling_number)+" out of "+str(total_reads)+"...")
    print("Found a total of "+str(total_reads)+" reads in blast file")

    spot_read=0
    zippy = 0
    for i in records:
        j=i["report"]["results"]["search"]
        try:
            spot_read+=1
            if spot_read%10000 == 0:
                print("Analyzing read: "+str(spot_read))
            hsps_number=0
            read_id = j["query_title"]
            read_length = j["query_len"]
            read_data_raw = j["hits"][0]["hsps"]
            for hsps in read_data_raw:
                read_data=hsps
                hsps_number+=1
                skipit="True"
                try:
                    if read_data['evalue'] < evalue_threshold:
                        skipit="False"
                except:
                    if "e-" in read_data["evalue"]:
                        skipit="False"
                try:
                    the_title = j["hits"][0]["description"][0]["title"]
                    if target_title == "":
                        skipit="False"
                    elif the_title.lower() == target_title.lower() and target_title != "":
                        skipit="False"
                        #print(the_title)
                    else:
                        skipit = "True"
                        #print(the_title)
                except:
                    continue

                h_from=read_data['hit_from']
                h_to=read_data['hit_to']
                q_from=read_data['query_from']
                q_to=read_data['query_to']
                q_strand=read_data['query_strand']
                h_strand=read_data['hit_strand']
                span=read_data['align_len']

                current_inserts = 0
                current_misses = 0
                current_dels = 0

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

                if skipit == "False":
                    hseq=read_data['hseq']
                    qseq=read_data['qseq']

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
                                    print(new_qseq[k].upper())
                                    print(genome_dict[str((ticker*k)+h_from)][2])
                                

                            #Tracks insertions and inserted bases
                            the_step = "insert"
                            if gymn_seq[0] == "-" and qseq[k] != "-":
                                genome_dict[str((ticker*k)+h_from)][5]+=1
                                try:
                                    genome_dict[str((ticker*k)+h_from)][6][qseq[k].upper()]+=1
                                except:
                                    print(qseq[k].upper())
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
                with open(save_locale+prefixed+'reads_df.csv', 'a') as f:
                    f.write(str(spot_read)+","+read_id+","+str(hsps_number)+","+str(read_data["bit_score"])+","+str(read_data["evalue"])+","+str(100*(span/len(genome)))+","+q_strand+","+h_strand+","+direction+","+str(read_length)+","+str(read_data['align_len'])+","+str(h_from)+","+str(h_to)+","+str(q_from)+","+str(q_to)+","+skipit+","+str(current_dels)+","+str(current_misses)+","+str(current_inserts)+"\n")
                max_insertion_length = 0
            if random_sample and aligned_reads >= sampling_number:
                break

        #If no alignment is found, the script passes
        except:
            unaligned_reads+=1
            pass
except:
    print("Yeah, something went wrong")

print("Aligned reads: "+str(aligned_reads)+" ("+str(100*round(aligned_reads/total_reads,2))+"%)")
if aligned_reads > 0:
    print("Aligned read errors: "+str(errors)+" ("+str(100*round(errors/aligned_reads,2))+"%)")
    print("\n")
    print("Generating genome df...")
    #print(genome_dict)
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
