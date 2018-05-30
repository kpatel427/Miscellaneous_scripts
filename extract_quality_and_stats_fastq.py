#!/usr/local/bin/python
# Khushbu Patel | 05/29/2018
#|__This script requires Python 3.4 and modules - numpy & scipy
#|__extracts the quality string and determine the length and average quality score of each read
#|__Converts the raw values for each read set into descriptive statistics
#|__Converts bam file to sam file format, extracts insert sizes of correctly paired mapped reads
#|__Provides descriptive stats for insert sizes 
# Usage: ./extract_quality_and_stats_fastq.py [input.fastq] [in.sam]

import numpy as np
import sys
from scipy.stats import skew,mstats


quality_scores = []
average_quality = 0
read_length = []
inserts = []
insert_sizes = []


fastq = sys.argv[1]
sam = sys.argv[2]

# To parse fastq file
def parseFastq(fastq_infile):
	sequences = []
	qualities = []
	
	with open(fastq_infile,"r", encoding="utf8", errors='ignore') as f:
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
	
	
# function call
seqs,quals = parseFastq(fastq)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals


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
	

	
# To parse input BAM file
def parseBam(input_bam):
	with open(input_bam,"r") as f:
		for line in f:
			if line.startswith('@'):	# skip the headers
				continue
			else:
				line = line.rstrip()
				temp = line.split()
				flag = temp[1]
				if(int(flag) == 99 or int(flag) == 163):
					insert_sizes.append(int(temp[8]))
		return insert_sizes
	
# function call
inserts = parseBam(sam)




# To calculate descriptive stats
def stats(in_array):
	a = np.array(in_array)
	mean = a.mean()
	std_dev = a.std()
	variance = np.var(a)
	Q1 = np.percentile(a,25)
	median = np.percentile(a,50)
	Q3 = np.percentile(a,75)
	skewness = skew(a)
	geometric_mean = mstats.gmean(a)
	
	high = []
	low = []
	IQR = Q3 - Q1
	lower = Q1 - (1.5*IQR)
	upper = Q3 - (1.5*IQR)
	
	if(min(in_array) < lower):
		low_whisker = min(in_array)
	else:
		low_whisker = min(in_array)
	
	if(max(in_array) > upper):
		high_whisker = max(in_array)
	else:
		high_whisker = max(in_array)
	
	return mean,std_dev,variance,Q1,median,Q3,skewness,geometric_mean,low_whisker,high_whisker
	
	
print("Number of reads for file is %d"  % len(seqs))	
	
	
	

# Descriptive stats for read length
r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read_length)
print("------ Descriptive stats for Read Lengths----------")
print("Normal Mean = %5.5f" % r_mean)
print("Geometric mean = %5.5f" % r_gmean)
print("standard deviation = %5.5f" % r_stdDev)
print("Variance = %5.5f" % r_var)
print("1st quartile = %5.5f" % r_Q1)
print("Median = %5.5f" % r_median)
print("3rd quartile = %5.5f" % r_Q3)
print("Skewness = %5.5f" % r_skew)
print("Lower whisker = %5.5f" % r_lwhisker)
print("upper whisker = %5.5f" % r_hwhisker)


# Descriptive stats for read quality
q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores)
print("------ Descriptive stats for Read Quality----------")
print("Normal Mean = %5.5f" % q_mean)
print("Geometric mean = %5.5f" % q_gmean)
print("standard deviation = %5.5f" % q_stdDev)
print("Variance = %5.5f" % q_var)
print("1st quartile = %5.5f" % q_Q1)
print("Median = %5.5f" % q_median)
print("3rd quartile = %5.5f" % q_Q3)
print("Skewness = %5.5f" % q_skew)
print("Lower whisker = %5.5f" % q_lwhisker)
print("upper whisker = %5.5f" % q_hwhisker)



# Descriptive stats for Insert sizes
i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(inserts)
print("------ Descriptive stats for Insert Sizes----------")
print("Normal Mean = %5.5f" % i_mean)
print("Geometric mean = %5.5f" % i_gmean)
print("standard deviation = %5.5f" % i_stdDev)
print("Variance = %5.5f" % i_var)
print("1st quartile = %5.5f" % i_Q1)
print("Median = %5.5f" % i_median)
print("3rd quartile = %5.5f" % i_Q3)
print("Skewness = %5.5f" % i_skew)
print("Lower whisker = %5.5f" % i_lwhisker)
print("upper whisker = %5.5f" % i_hwhisker)

