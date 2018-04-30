#!/usr/local/bin/python
# Khushbu Patel | 04/30/2018
# Calculates total number of bases in fasta file
# Usage: ./get_total_bases.py 

infile=""
infile=raw_input("Enter the file name:") #takes file as user input from the command line
f=open(infile,"r")
sum = 0

for line in f:
	line = line.rstrip()	#equivalent of chomp in perl; removes new line characters
	if line.startswith('>'):
		continue		#if starts with header, do nothing
	else:
		chars = list(line)	#parse the str to list
		sum += len(chars)	#count the list using len(); every iteration adds to the previous sum value
		
	
print(sum)	#print the sum
