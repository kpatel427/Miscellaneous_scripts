#!/usr/local/bin/python
# Khushbu Patel | 04/30/2018
# Calculates number of amibguous bases, percent ambiguous bases in fasta file
# Usage: ./get_ambi_bases.py 

infile=""
infile=raw_input("Enter the file name:") #takes file as user input from the command line
f=open(infile,"r")

ambi={}
ambi_codes = ['R','Y','S','W','K','M','B','D','H','V']

for line in f:
	line = line.rstrip()	#equivalent of chomp in perl; removes new line characters
	if line.startswith('>'):
		countN=0
		totalBases=0
		percentBases=0
		Id = line	#Assigning header to Id variable
		continue	
	else:
		countN += line.count("N")
		if any (x in ambi_codes for x in line):
			print ("There are ambiguous bases present, other than N")
		else:
			print ("No ambiguous bases other than N present!")
		ambi[Id] = countN
		chars = list(line)	#parse the string to list
		totalBases += len(chars)	#total number of bases in the contig
		percentBases = (countN/totalBases)*100
		
	
		
	
for x in ambi:		#printing the dictionary
	print(x,':',ambi[x],totalBases,percentBases)
	


