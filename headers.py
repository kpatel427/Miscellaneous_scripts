#!/usr/local/bin/python
# Khushbu Patel | 05/24/2018
# Usage: ./headers.py [input.fasta] [header.txt]

import sys

fasta_infile=sys.argv[1]
header_file = sys.argv[2]
new_fasta = "new_out.informative.fasta"

old_ID = []

x = open(new_fasta,"w")
f1 = open(header_file, "r")
f = open(fasta_infile,"r")

for line in f:
	if line.startswith('>'):
		line = line.rstrip()	# Removes trailing new lines
		temp = f1.readline()	# reads lines from file with new headers
		temp = temp.split()		# splits the line
		old_ID = temp[1]		# we want to replace old headers with information in column 2; storing each line into a list old_ID
		print >>x ,">"+old_ID	# writing the output to a new file
		
	else:
		line = line.rstrip()	# Removes trailing new lines
		print >> x, line		# writing the output to a new file

	

	
		
		

