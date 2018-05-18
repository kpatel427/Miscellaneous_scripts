#!/usr/local/bin/python
# Khushbu Patel | 05/15/2018
# Corrects the filter and Calculates number of bases falsely assigned incorrect filter when the bases pass the consensus and coverage. Prints out a new vcf with filters corrected. 
# Usage: ./correct_vcf_filter.py inputfile.vcf
# Output will be another vcf file, named correct_vcf.vcf; This script will not overwrite the original VCF

import sys

infile= sys.argv[1]	# Takes the file name as a command line argument
f=open(infile,"r")


temp = []
count =0
out = []
ID = ''

with open('correct_vcf.vcf', 'a') as f1:
	for line in f:
		line =line.rstrip()
		if line.startswith('#'):
			ID = line					# Printing the headers
			print >>f1, ID
		
		else:
			array = line.split()
			temp = array[9].split(':')
			
			if(array[6] != "PASS"):
				if(len(temp)> 6):
					temp[6] = temp[6].replace('%','')
					if(temp[3] >= 20 and float(temp[6]) < 5.0):		# if coverage and consensus both meet, change the filter to PASS
						#print(line, "--Incorrect Filter!")
						array[6] = "PASS"
						out = '\t'.join(array)
						count += 1					# Count the number of sites that have been assigned wrong filter
						print >>f1, out
					else:							# if allele frequency is greater than 5% or coverage does not meet; keep the filter
						out = '\t'.join(array)
						print >>f1, out

				else:									# else - print everything as is
					out = '\t'.join(array)
					print >>f1, out
				
			else:									# else - print everything as is
				out = '\t'.join(array)
				print >>f1, out
				
	print("Total sites lost due to incorrect filter %d" % count)	
		
