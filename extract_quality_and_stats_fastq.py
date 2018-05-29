#!/usr/local/bin/python
# Khushbu Patel | 05/29/2018
#extracts the quality string and determine the length and average quality score of each read; Converts the raw values for each read set into descriptive statistics
# Usage: ./extract_quality_and_stats_fastq.py [input.fastq]

import numpy as np
import sys

quality_scores = []
average_quality = 0
read_length = []


infile = sys.argv[1]

# To parse fastq file
def parseFastq(fastq_infile):
	sequences = []
	qualities = []
	
	with open(fastq_infile,"r") as f:
		while True:	
			f.readline()
			seq = f.readline().rstrip()		# gets sequence line
			f.readline()
			qual = f.readline().rstrip()	# gets quality line
			if len(seq) == 0:		# if seq length is 0; reached end of file so break out of the loop
				break	
			sequences.append(seq)	# append seq to sequences list
			qualities.append(qual)	# append qual to sequences list
	
	return sequences,qualities
	

seqs,quals = parseFastq(infile)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals


# To convert ASCII  to quality scores
def phred33toQ(qual):
	return ord(qual) - 33	# ord converts char to ASCII values and returns
	
	
# To get average quality scores for each read
for Q in quals:
	score = 0
	read_length.append(len(Q))
	for val in Q:
		score += phred33toQ(val)
	average_quality = (score/len(Q))	
	quality_scores.append(average_quality)		
	

	
# To calculate descriptive stats
def stats(in_array):
	a = np.array(in_array)
	mean = a.mean()
	std_dev = a.std()
	variance = np.var(a)
	Q1 = np.percentile(a,25)
	median = np.percentile(a,50)
	Q3 = np.percentile(a,75)
	
	return mean,std_dev,variance,Q1,median,Q3
	
print("Number of reads for file is %d"  % len(seqs))	
	
# Descriptive stats for read length
r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3 = stats(read_length)
print("------ Descriptive stats for Read Lengths----------")
print("Mean = %5.5f" % r_mean)
print("standard deviation = %5.5f" % r_stdDev)
print("Variance = %5.5f" % r_var)
print("1st quartile = %5.5f" % r_Q1)
print("Median = %5.5f" % r_median)
print("3rd quartile = %5.5f" % r_Q3)


# Descriptive stats for read quality
q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3 = stats(quality_scores)
print("------ Descriptive stats for Read Quality----------")
print("Mean = %5.5f" % q_mean)
print("standard deviation = %5.5f" % q_stdDev)
print("Variance = %5.5f" % q_var)
print("1st quartile = %5.5f" % q_Q1)
print("Median = %5.5f" % q_median)
print("3rd quartile = %5.5f" % q_Q3)

