import argparse
import uuid

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rnacanvas", help = "Path to the RNACanvas file.") #
parser.add_argument("-p", "--peFile", help = "Path to the PE file.") #
parser.add_argument("-d", "--divider", help = "Max float to divide PEs by.", default="2.7") #
parser.add_argument("-m", "--minOpacity", help = "Minimum opacity to be used.", default="0.05") #
parser.add_argument("-c", "--color", help = "Hex code for the fill color.", default='#ff0000') #
parser.add_argument("-o", "--output", help = "Name of the output file.") #
parser.parse_args()
args = parser.parse_args()
lcounter = 0
runIt=True

if args.output:
    outputFile=args.output
else:
    outputFile=args.rnacanvas.split(".rnacanvas")[0]+"_colored.rnacanvas"

try:
     minOpacity=float(args.minOpacity)
except:
    print("Given minimum opacity value unable to be read. Aborting")
    runIt=False

try:
     divider=float(args.divider)
except:
    print("Given divider unable to be read. Aborting")
    runIt=False

def load_data(peFile, canvasFile):
    with open(canvasFile, "r") as f:
        rnaCanvas = f.read()
    rnaCanvas.rstrip("\n")
    svg=rnaCanvas.split("</svg>")[0]
    front=svg.split("<svg")[0]
    back=rnaCanvas.split("</svg>")[1]
    mid=svg.split("<svg")[1]

    with open(peFile, "r") as f:
        pesAll = f.read()
    pesAll.strip("\n")
    pes=pesAll.split("\n")[2]
    seq=pesAll.split("\n")[1]
    pes=pes.replace("(", "")
    pes=pes.replace(")", "")
    tups=pes.split(" ")[1:]

    return front, mid, back, tups, seq

def stitch(front, mid, back):
    finalFile=front+"<svg"+mid+"</svg>"+back
    with open(outputFile, 'w') as f:
                f.write(finalFile)

def splitText(svg, pes, seq):
     gin={}
     starter=0
     text=svg.split("<text ")[1:]
     try:
        for i in range(0,len(seq)):
            liner=text[i].split("</text>")[0]
            x=liner.split(' x=')[1][2:]
            x=x.split('\"')[0]
            y=liner.split(' y=')[1][2:]
            y=y.split('\"')[0]
            gin[starter]=[liner[len(liner)-1], str(max(float(pes[starter])/divider, minOpacity)), str(float(x[:len(x)-1])+3), str(float(y[:len(y)-1])-3), uuid.uuid1()]
            starter+=1
     except:
          print("Unable to synch given PE values with given RNACanvas file.")
     return gin

def addCircles(text, svg, color):
     for i in text.keys():
          svg+="<circle id=\\\"uuid-"+str(text[i][4])+"\\\" cx=\\\""+text[i][2]+"\\\" cy=\\\""+text[i][3]+"\\\" stroke-width=\\\"0\\\" r=\\\"6\\\" fill=\\\""+color+"\\\" fill-opacity=\\\""+text[i][1]+"\\\"></circle>"
     return svg

if runIt:
    front, svg, back, pes, seq = load_data(args.peFile, args.rnacanvas)
    text = splitText(svg, pes, seq)
    if len(text.keys()) < len(seq):
         print("Skipping save step. Please double check your files.")
    else:
        svg = addCircles(text, svg, args.color)
        stitch(front, svg, back)