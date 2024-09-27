rna = "UUUAUGCUUAGUGAACUUGGCAUUAAGCAUGAA"
dotBra = ".((((((((((((.......))))))))))).)"
braDot = dotBra[::-1]

hitsLib = {'GC':0,
           'AU':0,
           'GU':0,
           'other':[]}
if len(dotBra) != len(rna):
    print("Sequences are not the same length. Aborting.")
elif dotBra.count('(') != dotBra.count(')'):
    print("Structure pair is not closed. Aborting.")
else:
    hits = 0
    for i in range(0,len(dotBra)):
        if dotBra[i] == "(":
            hits += 1
            stih = 0
            for j in range(0,len(braDot)):
                if braDot[j] == ')':
                    stih += 1

                if stih == hits:
                    pair = rna[i]+rna[::-1][j]
                    break

            print(pair)
            if pair in hitsLib.keys():
                hitsLib[pair]+=1
            elif pair[::-1] in hitsLib.keys():
                hitsLib[pair[::-1]]+=1
            else:
                hitsLib['other'].append(pair)

print(hitsLib)