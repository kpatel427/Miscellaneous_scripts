#!/usr/local/bin/python
# Khushbu Patel | 04/30/2018
# parses VCF
# Usage: ./parseVCF.py 

infile=""
infile=raw_input("Enter the file name:") #takes file as user input from the command line
f=open(infile,"r")

count =0
count_dp20 = 0
count_af = 0
count_pass = 0
count_other = 0


for line in f:
	line =line.rstrip()
	if line.startswith('##'):
		continue
	
	elif line.startswith('#CHROM'):
		header = line
		#print(line)
	
	else:
		array = line.split()
		#print(array)
		if(array[6] == 'DP20;AF0.95'):
			count = count+1
		elif(array[6] == 'DP20'):
			count_dp20 += 1
		elif(array[6] == 'AF0.95'):
			count_af += 1
		elif(array[6] == 'PASS'):
			count_pass += 1
		else:
			count_other +=1
			
			
print 'Number of bases that passed =', count_pass		
print 'Number of bases that did not meet both consensus and coverage =', count
print 'Number of bases that did not meet coverage only =', count_dp20
print 'Number of bases that did not meet consensus only =', count_af
print 'Number of bases that did not pass due to other reasons =', count_other
