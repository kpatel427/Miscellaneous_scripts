#!/usr/local/bin/python
# creating and counting occurences of kmers of length less than k from a given dna string of length s

l = "ACGTAG"
dna = "ACGGATCGGCATCGT"

count = {}
words = []


k = len(l)
s = len(dna)

for x in range(k-1):
  for i in range(0,s-x+1):
    if(dna[i:i+x] == ''):
      continue
    else:
      words.append(dna[i:i+x])



for a in words:
  if(a in count.keys()):
    count[a] = count[a] + 1
  else:
    count[a] = 1

for k,v in count.items():
  print(k,v)
