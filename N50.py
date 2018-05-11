#!/usr/local/bin/python
# Khushbu Patel | 05/11/2018
# Takes input a file with multiple contigs of varying lengths and gives N50 score (minimum contig length required to cover 50 percent of the assembled genome sequence)
# Usage: ./N50.py 


infile=""
infile=raw_input("Enter the file name:") #takes file as user input from the command line
f=open(infile,"r")

sum_contigs = 0
half = 0
n50 = []
lengths_contigs = []



for line in f:
	line = line.rstrip()
	sum_contigs += len(line)
	half =sum_contigs/2
	lengths_contigs.append(len(line))


n50 = sorted(lengths_contigs)	#sort contig lengths
len_n50 = len(n50)				#getting total number of elements in n50 list

if(len_n50 % 2 == 0):			#if even, avergage of two middle values is taken
	mid = (len_n50)/2
	avg = (mid+mid+1)/ 2
	print(n50[avg])
else:							#else, picks the middle value
	mid = (len_n50)/ 2
	print(n50[mid])
