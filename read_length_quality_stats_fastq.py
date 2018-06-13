#!/usr/local/bin/python
# Khushbu Patel | 05/29/2018
#|__This script requires Python 3.4 and modules - numpy & scipy
#|__extracts the quality string and determine the length and average quality score of each read
#|__Converts the raw values for each read set into descriptive statistics
#|__Provides descriptive stats for Read Lengths and Read Qualities, number and percentage of reads below Q30 and Ambiguous base counts
# Usage: ./extract_quality_and_stats_fastq.py [Read_1.fastq] [Read_2.fastq]

import numpy as np
import sys
from scipy.stats import skew,mstats

# ------------------------------------------ DECLARATIONS AND INITIALIZATIONS ------------------------------------------------#
quality_scores_R1 = []
quality_scores_R2 = []
average_quality = 0
read1_length = []
read2_length = []
inserts = []
insert_sizes = []
R1_le_249 = 0
R1_gt_249 = 0
R1_le_149 = 0
R1_gt_149 = 0
R1_le_299 = 0
R1_gt_299 = 0
R2_le_249 = 0
R2_gt_249 = 0
R2_le_149 = 0
R2_gt_149 = 0
R2_le_299 = 0
R2_gt_299 = 0
countN1 = 0
countN2 = 0
Q1_lt_30 = 0
Q2_lt_30 = 0
R1 = []
R2 = []
Q1 = []
Q2 = []

# ------------------------------------------ FUNCTIONS ------------------------------------------------#
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
	


# To convert ASCII  to quality scores
def phred33toQ(qual):
	return ord(qual) - 33	# ord converts char to ASCII values and returns
	


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
	
	
	
# Ambiguous base counts
def countN(seq):
	count = 0
	for s in seq:
		count += s.count("N")
	return count

	
# quality thresholds
def Q30(qual_list):
	count_lt_30 = 0
	for x in qual_list:
		if(x < 30):
			count_lt_30 += 1
		else:
			continue
	return count_lt_30
	
	
# To get average quality scores for each read1 
def qual_score(qual):
	quality_scores = []
	read_len = []
	for Q in qual:
		score = 0
		read_len.append(len(Q))
		for val in Q:
			score += phred33toQ(val)
		average_quality = (score/len(Q))	
		quality_scores.append(average_quality)	
	return read_len,quality_scores
	

# read Length thresholds 
def threshold(lengths):
	len_le_149 = 0
	len_gt_149 = 0
	len_le_249 = 0
	len_gt_249 = 0
	len_le_299 = 0
	len_gt_299 = 0
	for x in lengths:
		if(x <= 149):
			len_le_149 += 1
		elif(x > 149):
			len_gt_149 += 1
		elif(x <= 249):
			len_le_249 += 1
		elif(x > 249):
			len_gt_249 += 1
		elif(x <= 299):
			len_le_299 += 1
		elif(x > 299):
			len_gt_299 += 1

	return len_le_149,len_gt_149,len_le_249,len_gt_249,len_le_299,len_gt_299
	
	
	
# ---------------------------------------------------- MAIN ----------------------------------------------------------------- #	

# command line arguments
fastq1 = sys.argv[1]
fastq2 = sys.argv[2]

	

# Parsing fastq: function call
seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
seqs2,quals2 = parseFastq(fastq2)
	
# total number of reads
read_count1 = len(seqs1)
read_count2 = len(seqs2)


# average quality scores for each read: function call
read1_length,quality_scores_R1 = qual_score(quals1)
read2_length,quality_scores_R2 = qual_score(quals2)


# read Length thresholds: function call
R1_le_149,R1_gt_149,R1_le_249,R1_gt_249,R1_le_299,R1_gt_299 = threshold(read1_length)
R2_le_149,R2_gt_149,R2_le_249,R2_gt_249,R2_le_299,R2_gt_299 = threshold(read2_length)

# quality threshold function call: function call
Q1_lt_30 = Q30(quality_scores_R1)
Q2_lt_30 = Q30(quality_scores_R2)

percent_reads_lt_30_R1 = Q1_lt_30/len(seqs1) * 100
percent_reads_lt_30_R2 = Q2_lt_30/len(seqs2) * 100

# Ambiguous base function call: function call
countN1 = countN(seqs1)
countN2 = countN(seqs2)

# Descriptive stats for read1 length: function call
r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read1_length)
i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(read2_length)


# Descriptive stats for Q1 quality: function call
q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores_R1)
s_mean,s_stdDev,s_var,s_Q1,s_median,s_Q3,s_skew,s_gmean,s_lwhisker,s_hwhisker = stats(quality_scores_R2)


# Result lists
if(R1_le_149 > 0 and R1_gt_149 > 0):
	perc_R1_le_149 = (R1_le_149/read_count1) * 100
	perc_R1_gt_149 = (R1_gt_149/read_count1) * 100
	perc_R2_le_149 = (R2_le_149/read_count2) * 100
	perc_R2_gt_149 = (R2_gt_149/read_count2) * 100
	
	R1 = [["\t",'Stats for R1 Length',"\t",'Stats for R2 Length'],['Normal mean:	',r_mean,"\t",i_mean],['SD:	',r_stdDev,"\t",i_stdDev],['Variance:	',r_var,"\t",i_var],['1st quartile:	',r_Q1,"\t",i_Q1],['Median:	',r_median,"\t",i_median],
			['3rd quartile:	',r_Q3,"\t",i_Q3],['Skewness:	',r_skew,"\t",i_skew],['Geormetric mean:	',r_gmean,"\t",i_gmean],['Lower whisker:	',r_lwhisker,"\t",i_lwhisker],['Higher whisker:	',r_hwhisker,"\t",i_hwhisker],
			['Reads <= 149:	',R1_le_149,"\t",R2_le_149],['Reads > 149:	',R1_gt_149,"\t",R2_gt_149],['Percent counts for reads <= 149:	',perc_R1_le_149,"\t",perc_R2_le_149],['Percent counts for reads > 149:	',perc_R1_gt_149,"\t",perc_R2_gt_149]]


	Q1 = [["\t",'Stats for R1 Quality',"\t",'Stats for R2 Quality'],['Normal mean:	',q_mean,"\t",s_mean],['SD:	',q_stdDev,"\t",s_stdDev],['Variance:	',q_var,"\t",s_var],['1st quartile:	',q_Q1,"\t",s_Q1],['Median:	',q_median,"\t",s_median],['3rd quartile:	',q_Q3,"\t",s_Q3],
			['Skewness:	',q_skew,"\t",s_skew],['Geometric mean:	',q_gmean,"\t",s_gmean],['Lower whisker:	',q_lwhisker,"\t",s_lwhisker],
			['Higher whisker:	',q_hwhisker,"\t",s_hwhisker],['Reads count:	',read_count1,"\t",read_count2],['Reads below Q30:	',Q1_lt_30,"\t",Q2_lt_30],['Percentage of reads below Q30:	',percent_reads_lt_30_R1,"\t",percent_reads_lt_30_R2],
			['Ambiguous bases:	',countN1,"\t",countN2]]
			
elif(R1_le_249 > 0 and R1_gt_249 > 0 ):
	perc_R1_le_249 = (R1_le_249/read_count1) * 100
	perc_R1_gt_249 = (R1_gt_249/read_count1) * 100
	perc_R2_le_249 = (R2_le_249/read_count2) * 100
	perc_R2_gt_249 = (R2_gt_249/read_count2) * 100
	
	R1 = [["\t",'Stats for R1 Length',"\t",'Stats for R2 Length'],['Normal mean:	',r_mean,"\t",i_mean],['SD:	',r_stdDev,"\t",i_stdDev],['Variance:	',r_var,"\t",i_var],['1st quartile:	',r_Q1,"\t",i_Q1],['Median:	',r_median,"\t",i_median],
			['3rd quartile:	',r_Q3,"\t",i_Q3],['Skewness:	',r_skew,"\t",i_skew],['Geormetric mean:	',r_gmean,"\t",i_gmean],['Lower whisker:	',r_lwhisker,"\t",i_lwhisker],['Higher whisker:	',r_hwhisker,"\t",i_hwhisker],
			['Reads <= 249:	',R1_le_249,"\t",R2_le_249],['Reads > 249:	',R1_gt_249,"\t",R2_gt_249],['Percent counts for reads <= 249:	',perc_R1_le_249,"\t",perc_R2_le_249],['Percent counts for reads > 249:	',perc_R1_gt_249,"\t",perc_R2_gt_249]]


	Q1 = [["\t",'Stats for R1 Quality',"\t",'Stats for R2 Quality'],['Normal mean:	',q_mean,"\t",s_mean],['SD:	',q_stdDev,"\t",s_stdDev],['Variance:	',q_var,"\t",s_var],['1st quartile:	',q_Q1,"\t",s_Q1],['Median:	',q_median,"\t",s_median],['3rd quartile:	',q_Q3,"\t",s_Q3],
			['Skewness:	',q_skew,"\t",s_skew],['Geometric mean:	',q_gmean,"\t",s_gmean],['Lower whisker:	',q_lwhisker,"\t",s_lwhisker],
			['Higher whisker:	',q_hwhisker,"\t",s_hwhisker],['Reads count:	',read_count1,"\t",read_count2],['Reads below Q30:	',Q1_lt_30,"\t",Q2_lt_30],['Percentage of reads below Q30:	',percent_reads_lt_30_R1,"\t",percent_reads_lt_30_R2],
			['Ambiguous bases:	',countN1,"\t",countN2]]
			
else:
	perc_R1_le_299 = (R1_le_299/read_count1) * 100
	perc_R1_gt_299 = (R1_gt_299/read_count1) * 100
	perc_R2_le_299 = (R2_le_299/read_count2) * 100
	perc_R2_gt_299 = (R2_gt_299/read_count2) * 100
	
	R1 = [["\t",'Stats for R1 Length',"\t",'Stats for R2 Length'],['Normal mean:	',r_mean,"\t",i_mean],['SD:	',r_stdDev,"\t",i_stdDev],['Variance:	',r_var,"\t",i_var],['1st quartile:	',r_Q1,"\t",i_Q1],['Median:	',r_median,"\t",i_median],
			['3rd quartile:	',r_Q3,"\t",i_Q3],['Skewness:	',r_skew,"\t",i_skew],['Geormetric mean:	',r_gmean,"\t",i_gmean],['Lower whisker:	',r_lwhisker,"\t",i_lwhisker],['Higher whisker:	',r_hwhisker,"\t",i_hwhisker],
			['Reads <= 299:	',R1_le_299,"\t",R2_le_299],['Reads > 299:	',R1_gt_299,"\t",R2_gt_299],['Percent counts for reads <= 299:	',perc_R1_le_299,"\t",perc_R2_le_299],['Percent counts for reads > 299:	',perc_R1_gt_299,"\t",perc_R2_gt_299]]


	Q1 = [["\t",'Stats for R1 Quality',"\t",'Stats for R2 Quality'],['Normal mean:	',q_mean,"\t",s_mean],['SD:	',q_stdDev,"\t",s_stdDev],['Variance:	',q_var,"\t",s_var],['1st quartile:	',q_Q1,"\t",s_Q1],['Median:	',q_median,"\t",s_median],['3rd quartile:	',q_Q3,"\t",s_Q3],
			['Skewness:	',q_skew,"\t",s_skew],['Geometric mean:	',q_gmean,"\t",s_gmean],['Lower whisker:	',q_lwhisker,"\t",s_lwhisker],
			['Higher whisker:	',q_hwhisker,"\t",s_hwhisker],['Reads count:	',read_count1,"\t",read_count2],['Reads below Q30:	',Q1_lt_30,"\t",Q2_lt_30],['Percentage of reads below Q30:	',percent_reads_lt_30_R1,"\t",percent_reads_lt_30_R2],
			['Ambiguous bases:	',countN1,"\t",countN2]]

	
for a,b,c,d in R1:
	print('%s %s %s %s' %(a, b, c, d))

print("\n-----------------------------------------------------------\n")
	
for a,b,c,d in Q1:
	print('%s %s %s %s' %(a, b, c, d))	

