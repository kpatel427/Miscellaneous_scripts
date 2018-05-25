#!/usr/local/bin/python
# Khushbu Patel | 05/25/2018
# Usage: ./parseFastq.py reads.fastq

import sys

sequences = []
qualities = []
fastq_infile=sys.argv[1]

with open(fastq_infile,"r") as f:
	while True:	
		f.readline()
		seq = f.readline().rstrip()
		f.readline()
		qual = f.readline().rstrip()
		if len(seq) == 0:
			break
		sequences.append(seq)
		qualities.append(qual)


print(sequences)
print(qualities)
