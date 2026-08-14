import argparse
import random
import re
import RNA

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--structure", help = "Structure to be copied", default="(((((((...)))))))") #
parser.add_argument('-p', '--percentGC', help="Percent (0.0-1.0) GC of the random bases.", default=0.5)
parser.add_argument('-c', '--coding', help= "Keep coding sequence the same. Should be an integer indicating the position where coding begins. Default is 0 (first base)")
parser.add_argument('-e', '--endCoding', help= "Position to end the coding on.")
parser.add_argument('-u', '--usageTable', help = 'Path to codon usage table')
parser.add_argument('-F', '--Force', help = 'Whether new sequence should be folded to verify structure. Default is False', default=False, action='store_true')
parser.add_argument('-n', '--numberOfIterations', help="Number of iterations to make.", default=100)
parser.add_argument("-o", "--out", help="Name of output file", default="mimicked.fasta") #
parser.parse_args()
args = parser.parse_args()
runIt = True

codon_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
        }

pairings = {
        "A":"U",
        "U":"A",
        "C":"G",
        "G":"C"
        }
def pickNucleotide(percentGC):
        position=random.randint(1,100)/100
        if position <= float(percentGC):
                return random.choice(['G','C'])
        else:
                return random.choice(['U', 'A'])

apical="C[U]+G"
targetStructure=args.structure.replace(".", "U")
targetStructure=targetStructure.replace("(", "C")
targetStructure=targetStructure.replace(")", "G")

firstHalf=re.split(apical, targetStructure)[0]+"C"
secondHalf="G"+re.split(apical, targetStructure)[1]
apical=targetStructure[len(firstHalf):(len(targetStructure)-len(secondHalf))]

if runIt:
        attempts=0
        for n in range(0,int(args.numberOfIterations)):
                firstHalfTemp=firstHalf

                #Randomly generate a 5' sequence
                newFirst=""
                for i in firstHalf:
                        newFirst+=pickNucleotide(float(args.percentGC))
                #Randomly generate the apical loop sequence
                newApical=""
                for i in range(0, len(apical)):
                        if i < len(apical)-1:
                                newApical+=pickNucleotide(float(args.percentGC))
                        else:
                                trick=pickNucleotide(float(args.percentGC))
                                while pairings[trick] != newApical[0]:
                                        trick=pickNucleotide(float(args.percentGC))
                                newApical+=trick
                newSeq=newFirst+newApical
                #Fill in the 3' sequence based on the 5' sequence/structure
                newSecond=""
                for i in secondHalf:
                        if len(firstHalfTemp) > 0:
                                while firstHalfTemp[-1] == "U":
                                        firstHalfTemp=firstHalfTemp[:len(firstHalfTemp)-1]
                                        newFirst=newFirst[:len(newFirst)-1]
                                if i =="U" and firstHalfTemp[-1] != "U":      
                                        newSecond+=pickNucleotide(float(args.percentGC))
                                elif (i == "C" and firstHalfTemp[-1] == "G") or (i == "G" and firstHalfTemp[-1] == "C"):
                                        newSecond+=pairings[newFirst[-1]]
                                        firstHalfTemp=firstHalfTemp[:len(firstHalfTemp)-1]
                                        newFirst=newFirst[:len(newFirst)-1]
                        else:
                                newSecond+=pickNucleotide(float(args.percentGC))
                if args.Force:
                        if RNA.fold(newSeq+newSecond)[0] != args.structure:
                                #print("Structure not maintained. Skipping sequence.")
                                attempts+=1
                                if attempts%(10*int(args.numberOfIterations)) == 0:  # Limit the number of failed attempts
                                        print("Failed "+str(attempts)+" times. Stopping.")      
                                        break
                                continue
                outText=">Mimicked Sequence_"+str(n+1)+"\n"
                outText+=newSeq+newSecond+"\n"
                with open(args.out, 'a') as f:
                        f.write(outText)
