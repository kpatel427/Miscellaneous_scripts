#!/usr/local/bin/python
# Khushbu Patel | 04/30/2018
# Calculates number of amibguous bases, percent ambiguous bases in fasta file
# Usage: ./get_ambi_bases.py 

infile=""
infile=raw_input("Enter the file name:") #takes file as user input from the command line
f=open(infile,"r")


ambi_codes = ['R','Y','S','W','K','M','B','D','H','V']
countN=0
Id = ''
totalBases=0
percentBases=0
other = 0

for line in f:
	line = line.rstrip()	#equivalent of chomp in perl; removes new line characters
	if line.startswith('>'):
		Id = line	#Assigning header to Id variable
		continue	
	else:
		countN += line.count("N")
		if any (x in ambi_codes for x in line):
			other += 1
		else:
			continue
		totalBases += len(line)	#total number of bases 
		percentBases = (countN/totalBases)*100
		
	
		
	
print ("Total number of N's = %d" % countN)
print ("Percentage of N's = %d" % percentBases)
print ("Total number of ambi bases other than N's = %d" % other)

	


