#!/usr/local/bin/python
# Khushbu Patel | 05/29/2018
#|__This script requires Python 3.4 and modules - numpy & scipy
#|__extracts the quality string and determine the length and average quality score of each read
#|__Converts the raw values for each read set into descriptive statistics
#|__Provides descriptive stats for Read Lengths and Read Qualities, number and percentage of reads below Q30 and Ambiguous base counts
#|__Outputs separate tables for different read length buckets (150bp,250bp and 300bp)
# Usage: ./read_length_quality_and_stats_fastq.py 

import numpy as np
from scipy.stats import skew,mstats
import glob


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
dic_149_bucket = {}
dic_249_bucket = {}
dic_299_bucket = {}
file1 = []
file2 = []
tot_len1_gt_149 = 0
tot_len1_le_249 = 0
tot_len1_gt_249 = 0
tot_len1_le_299 = 0
tot_len1_gt_299 = 0
tot_len2_gt_149 = 0
tot_len2_le_249 = 0
tot_len2_gt_249 = 0
tot_len2_le_299 = 0
tot_len2_gt_299 = 0



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
	tot_len_le_149 = 0
	tot_len_gt_149 = 0
	tot_len_le_249 = 0
	tot_len_gt_249 = 0
	tot_len_le_299 = 0
	tot_len_gt_299 = 0
	
	for x in lengths:
		if(x <= 149):
			len_le_149 += 1		# number of reads less than 149
			tot_len_le_149 += x	# Total length of all reads less than 149
		elif(x < 245):
			len_gt_149 += 1
			tot_len_gt_149 += x
		elif(x <= 249):
			len_le_249 += 1
			tot_len_le_249 += x
		elif(x < 295):
			len_gt_249 += 1
			tot_len_gt_249 += x
		elif(x <= 299):
			len_le_299 += 1
			tot_len_le_299 += x
		elif(x > 299):
			len_gt_299 += 1
			tot_len_gt_299 += x

	return len_le_149,len_gt_149,len_le_249,len_gt_249,len_le_299,len_gt_299,tot_len_le_149,tot_len_gt_149,tot_len_le_249,tot_len_gt_249,tot_len_le_299,tot_len_gt_299
	



def avg_qual_score(read_length,quality_score): 

	for l,q in zip(read_length,quality_score):

		if(l >= 145 and l < 245): # for lengths above and below 149
			dic_149_bucket[l] = q
		elif(l >= 245 and l < 295): # for lengths above and below 249	
			dic_249_bucket[l] = q
		elif(l >= 295):	# for lengths equal to or greater than 299
			dic_299_bucket[l] = q
		

	sum = 0
	num = 0
	avg_quality = 0

	
	if(len(dic_149_bucket) > 0):
		for k,v in dic_149_bucket.items():
			sum += v
			num += 1
			avg_quality = (sum/num)
		
		return avg_quality
		
	elif(len(dic_249_bucket) > 0):
		for k,v in dic_249_bucket.items():
			sum += v
			num += 1
			avg_quality = (sum/num)
			
		return avg_quality
		
	elif(len(dic_299_bucket) > 0):
		for k,v in dic_299_bucket.items():
			sum += v
			num += 1
			avg_quality = (sum/num)
			
		return avg_quality
	

	
	
def print_150bp():	
	print("\n\n-----Stats for 149 bucket---------")
	print(*files_149, sep='\t')
	print("Read Length Stats:")
	print(*N_mean_149, sep='\t')
	print(*G_mean_149, sep='\t')
	print(*SD_149, sep='\t')
	print(*Variance_149, sep='\t')
	print(*Skew_149, sep='\t')
	print(*median_149, sep='\t')
	print(*Q1_149, sep='\t')
	print(*Q3_149, sep='\t')
	print(*lwhisker_149, sep='\t')
	print(*hwhisker_149, sep='\t')
	print(*R_le_149, sep='\t')
	print(*R_gt_149, sep='\t')
	print(*final_perc_R1_le_149, sep='\t')
	print(*final_perc_R1_gt_149, sep='\t')
	print(*final_avg_quality_149, sep='\t')
	print(*final_avg_length_149, sep='\t')
	print("\nRead Quality Stats:")
	print(*qual_N_mean_149, sep='\t')
	print(*qual_G_mean_149, sep='\t')
	print(*qual_SD_149, sep='\t')
	print(*qual_Variance_149, sep='\t')
	print(*qual_skew_149, sep='\t')
	print(*qual_median_149, sep='\t')
	print(*qual_Q1_149, sep='\t')
	print(*qual_Q3_149, sep='\t')
	print(*qual_lwhisker_149, sep='\t')
	print(*qual_hwhisker_149, sep='\t')
	print(*total_no_reads_149, sep='\t')
	print(*qual_lt_30_149, sep='\t')
	print(*perc_qual_lt_30_149, sep='\t')
	print(*ambi_calls_149, sep='\t')
	
def print_250bp():	
	print("\n\n-----Stats for 249 bucket---------")
	print(*files_249, sep='\t')
	print("Read Length Stats:")
	print(*N_mean_249, sep='\t')
	print(*G_mean_249, sep='\t')
	print(*SD_249, sep='\t')
	print(*Variance_249, sep='\t')
	print(*Skew_249, sep='\t')
	print(*median_249, sep='\t')
	print(*Q1_249, sep='\t')
	print(*Q3_249, sep='\t')
	print(*lwhisker_249, sep='\t')
	print(*hwhisker_249, sep='\t')
	print(*R_le_249, sep='\t')
	print(*R_gt_249, sep='\t')
	print(*final_perc_R1_le_249, sep='\t')
	print(*final_perc_R1_gt_249, sep='\t')
	print(*final_avg_quality_249, sep='\t')
	print(*final_avg_length_249, sep='\t')
	print("\nRead Quality Stats:")
	print(*qual_N_mean_249, sep='\t')
	print(*qual_G_mean_249, sep='\t')
	print(*qual_SD_249, sep='\t')
	print(*qual_Variance_249, sep='\t')
	print(*qual_skew_249, sep='\t')
	print(*qual_median_249, sep='\t')
	print(*qual_Q1_249, sep='\t')
	print(*qual_Q3_249, sep='\t')
	print(*qual_lwhisker_249, sep='\t')
	print(*qual_hwhisker_249, sep='\t')
	print(*total_no_reads_249, sep='\t')
	print(*qual_lt_30_249, sep='\t')
	print(*perc_qual_lt_30_249, sep='\t')
	print(*ambi_calls_249, sep='\t')

	
def print_300bp():	
	print("\n\n-----Stats for 299 bucket---------")
	print(*files_299, sep='\t')
	print("Read Length Stats:")
	print(*N_mean_299, sep='\t')
	print(*G_mean_299, sep='\t')
	print(*SD_299, sep='\t')
	print(*Variance_299, sep='\t')
	print(*Skew_299, sep='\t')
	print(*median_299, sep='\t')
	print(*Q1_299, sep='\t')
	print(*Q3_299, sep='\t')
	print(*lwhisker_299, sep='\t')
	print(*hwhisker_299, sep='\t')
	print(*R_le_299, sep='\t')
	print(*R_gt_299, sep='\t')
	print(*final_perc_R1_le_299, sep='\t')
	print(*final_perc_R1_gt_299, sep='\t')
	print(*final_avg_quality_299, sep='\t')
	print(*final_avg_length_299, sep='\t')
	print("\nRead Quality Stats:")
	print(*qual_N_mean_299, sep='\t')
	print(*qual_G_mean_299, sep='\t')
	print(*qual_SD_299, sep='\t')
	print(*qual_Variance_299, sep='\t')
	print(*qual_skew_299, sep='\t')
	print(*qual_median_299, sep='\t')
	print(*qual_Q1_299, sep='\t')
	print(*qual_Q3_299, sep='\t')
	print(*qual_lwhisker_299, sep='\t')
	print(*qual_hwhisker_299, sep='\t')
	print(*total_no_reads_299, sep='\t')
	print(*qual_lt_30_299, sep='\t')
	print(*perc_qual_lt_30_299, sep='\t')
	print(*ambi_calls_299, sep='\t')

			
	
# ---------------------------------------------------- MAIN ----------------------------------------------------------------- #	

files_149 = [] #Stores paired read files
files_249 = [] #Stores paired read files
files_299 = [] #Stores paired read files
# Following lists are to store all results for 149bp bucket
N_mean_149 = ["Mean:"]
SD_149 = ["Std Deviation:"]
Variance_149 = ["Variance"]
median_149 = ["Median"]
Q1_149 = ["1st Quartile:"]
Q3_149 = ["3rd Quartile:"]
lwhisker_149 = ["Lower whisker:"]
hwhisker_149 = ["Upper Whisker:"]
Skew_149 = ["Skewness:"]
G_mean_149 = ["Geometric Mean:"]

qual_N_mean_149 = ["Mean:"]
qual_SD_149 = ["Std Deviation:"]
qual_Variance_149 = ["Variance:"]
qual_median_149 = ["Median:"]
qual_Q1_149 = ["1st Quartile:"]
qual_Q3_149 = ["3rd Quartile:"]
qual_lwhisker_149 = ["Lower whisker:"]
qual_hwhisker_149 = ["Upper Whisker:"]
qual_skew_149 = ["Skewness:"]
qual_G_mean_149 = ["Geometric Mean:"]

# Following lists are to store all results for 249bp bucket
N_mean_249 = ["Mean:"]
SD_249 = ["Std Deviation:"]
Variance_249 = ["Variance"]
median_249 = ["Median"]
Q1_249 = ["1st Quartile:"]
Q3_249 = ["3rd Quartile:"]
lwhisker_249 = ["Lower whisker:"]
hwhisker_249 = ["Upper Whisker:"]
Skew_249 = ["Skewness:"]
G_mean_249 = ["Geometric Mean:"]

qual_N_mean_249 = ["Mean:"]
qual_SD_249 = ["Std Deviation:"]
qual_Variance_249 = ["Variance:"]
qual_median_249 = ["Median:"]
qual_Q1_249 = ["1st Quartile:"]
qual_Q3_249 = ["3rd Quartile:"]
qual_lwhisker_249 = ["Lower whisker:"]
qual_hwhisker_249 = ["Upper Whisker:"]
qual_skew_249 = ["Skewness:"]
qual_G_mean_249 = ["Geometric Mean:"]

# Following lists are to store all results for 299bp bucket
N_mean_299 = ["Mean:"]
SD_299 = ["Std Deviation:"]
Variance_299 = ["Variance"]
median_299 = ["Median"]
Q1_299 = ["1st Quartile:"]
Q3_299 = ["3rd Quartile:"]
lwhisker_299 = ["Lower whisker:"]
hwhisker_299 = ["Upper Whisker:"]
Skew_299 = ["Skewness:"]
G_mean_299 = ["Geometric Mean:"]

qual_N_mean_299 = ["Mean:"]
qual_SD_299 = ["Std Deviation:"]
qual_Variance_299 = ["Variance:"]
qual_median_299 = ["Median:"]
qual_Q1_299 = ["1st Quartile:"]
qual_Q3_299 = ["3rd Quartile:"]
qual_lwhisker_299 = ["Lower whisker:"]
qual_hwhisker_299 = ["Upper Whisker:"]
qual_skew_299 = ["Skewness:"]
qual_G_mean_299 = ["Geometric Mean:"]

total_no_reads_149 = ["Read count:"]
total_no_reads_249 = ["Read count:"]
total_no_reads_299 = ["Read count:"]
qual_lt_30_149 = ["Reads < Q30:"]
qual_lt_30_249 = ["Reads < Q30:"]
qual_lt_30_299 = ["Reads < Q30:"]
perc_qual_lt_30_149 = ["Percentage reads < Q30"]
perc_qual_lt_30_249 = ["Percentage reads < Q30"]
perc_qual_lt_30_299 = ["Percentage reads < Q30"]
ambi_calls_149 = ["Amibguous base calls:"]
ambi_calls_249 = ["Amibguous base calls:"]
ambi_calls_299 = ["Amibguous base calls:"]


R_le_149 = ["Reads <= 149:"]
R_gt_149 = ["Reads > 149:"]
R_le_249 = ["Reads <= 249:"]
R_gt_249 = ["Reads > 249:"]
R_le_299 = ["Reads <= 299:"]
R_gt_299 = ["Reads > 299:"]

r_median = 0
i_median = 0

final_perc_R1_le_149 = ["% Reads <= 149:"]
final_perc_R1_gt_149 = ["% Reads > 149:"]


final_perc_R1_le_249 = ["% Reads <= 249:"]
final_perc_R1_gt_249 = ["% Reads > 249:"]


final_perc_R1_le_299 = ["% Reads <= 299:"]
final_perc_R1_gt_299 = ["% Reads > 299:"]


final_avg_quality_149 = ["Avergage Quality:"]
final_avg_length_149 = ["Average Length"]
final_avg_quality_249 = ["Avergage Quality:"]
final_avg_length_249 = ["Average Length"]
final_avg_quality_299 = ["Avergage Quality:"]
final_avg_length_299 = ["Average Length"]




# getting all fastq files
temp1 = []

for name in glob.glob('./*_1.fastq'):
    file1.append(name)

for name in glob.glob('./*_2.fastq'):
    file2.append(name)
	
file1 = sorted(file1)
file2 = sorted(file2)

for f1,f2 in zip(file1,file2):
	#print(f1,f2)
	
	#files.extend((f1,f2))
	

	# Result lists
	if(r_median < 152 and i_median < 152):
		files_149.extend((f1,f2))
		
		# command line arguments
		fastq1 = f1
		fastq2 = f2

		# Parsing fastq: function call
		seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
		seqs2,quals2 = parseFastq(fastq2)
			
		# total number of reads
		read_count1 = len(seqs1)
		read_count2 = len(seqs2)
		
		total_no_reads_149.extend((read_count1,read_count2))


		# average quality scores for each read: function call
		read1_length,quality_scores_R1 = qual_score(quals1)
		read2_length,quality_scores_R2 = qual_score(quals2)
		
		
		# quality threshold function call: function call
		Q1_lt_30 = Q30(quality_scores_R1)
		Q2_lt_30 = Q30(quality_scores_R2)
		qual_lt_30_149.extend((Q1_lt_30,Q2_lt_30))

		percent_reads_lt_30_R1 = Q1_lt_30/len(seqs1) * 100
		percent_reads_lt_30_R2 = Q2_lt_30/len(seqs2) * 100
		perc_qual_lt_30_149.extend((percent_reads_lt_30_R1,percent_reads_lt_30_R2))

		# Ambiguous base function call: function call
		countN1 = countN(seqs1)
		countN2 = countN(seqs2)
		ambi_calls_149.extend((countN1,countN2))
		
		# Descriptive stats for read1 length: function call
		r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read1_length)
		i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(read2_length)
		
		N_mean_149.extend((r_mean,i_mean))
		SD_149.extend((r_stdDev,i_stdDev))
		Variance_149.extend((r_var,i_var))
		median_149.extend((r_median,i_median))
		Q1_149.extend((r_Q1,i_Q1))
		Q3_149.extend((r_Q3,i_Q3))
		lwhisker_149.extend((r_lwhisker,i_lwhisker))
		hwhisker_149.extend((r_hwhisker,i_hwhisker))
		Skew_149.extend((r_skew,i_skew))
		G_mean_149.extend((r_gmean,i_gmean))

		
		# Descriptive stats for Q1 quality: function call
		q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores_R1)
		s_mean,s_stdDev,s_var,s_Q1,s_median,s_Q3,s_skew,s_gmean,s_lwhisker,s_hwhisker = stats(quality_scores_R2)
		
		qual_N_mean_149.extend((q_mean,s_mean))
		qual_SD_149.extend((q_stdDev,s_stdDev))
		qual_Variance_149.extend((q_var,s_var))
		qual_median_149.extend((q_median,s_median))
		qual_Q1_149.extend((q_Q1,s_Q1))
		qual_Q3_149.extend((q_Q3,s_Q3))
		qual_lwhisker_149.extend((q_lwhisker,s_lwhisker))
		qual_hwhisker_149.extend((q_hwhisker,s_hwhisker))
		qual_skew_149.extend((q_skew,s_skew))
		qual_G_mean_149.extend((q_gmean,s_gmean))
		
		
		

		# read Length thresholds: function call
		R1_le_149,R1_gt_149,R1_le_249,R1_gt_249,R1_le_299,R1_gt_299,tot_len1_le_149,tot_len1_gt_149,tot_len1_le_249,tot_len1_gt_249,tot_len1_le_299,tot_len1_gt_299 = threshold(read1_length)
		R2_le_149,R2_gt_149,R2_le_249,R2_gt_249,R2_le_299,R2_gt_299,tot_len2_le_149,tot_len2_gt_149,tot_len2_le_249,tot_len2_gt_249,tot_len2_le_299,tot_len2_gt_299 = threshold(read2_length)
		
		
		R_le_149.extend((R1_le_149,R2_le_149))
		R_gt_149.extend((R1_gt_149,R2_gt_149))
		
		
		perc_R1_le_149 = (R1_le_149/read_count1) * 100
		perc_R1_gt_149 = (R1_gt_149/read_count1) * 100
		perc_R2_le_149 = (R2_le_149/read_count2) * 100
		perc_R2_gt_149 = (R2_gt_149/read_count2) * 100
	
		final_perc_R1_le_149.extend((perc_R1_le_149,perc_R2_le_149))
		final_perc_R1_gt_149.extend((perc_R1_gt_149,perc_R2_gt_149))
		
	
		# avg_qual_score: function call
		avg_quality_1 = avg_qual_score(read1_length,quality_scores_R1)
		avg_quality_2 = avg_qual_score(read2_length,quality_scores_R2)
		
		final_avg_quality_149.extend((avg_quality_1,avg_quality_2))
		
		avg_length_1 = (tot_len1_le_149 + tot_len1_gt_149) / (R1_le_149 + R1_gt_149)
		avg_length_2 = (tot_len2_le_149 + tot_len2_gt_149) / (R2_le_149 + R2_gt_149)
		
		final_avg_length_149.extend((avg_length_1,avg_length_2))
		
		
				
	elif(r_median < 252 and i_median < 252 ):
		files_249.extend((f1,f2))
		
		
		# command line arguments
		fastq1 = f1
		fastq2 = f2

		# Parsing fastq: function call
		seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
		seqs2,quals2 = parseFastq(fastq2)
			
		# total number of reads
		read_count1 = len(seqs1)
		read_count2 = len(seqs2)
		
		total_no_reads_249.extend((read_count1,read_count2))
		
		# average quality scores for each read: function call
		read1_length,quality_scores_R1 = qual_score(quals1)
		read2_length,quality_scores_R2 = qual_score(quals2)
		
		
		# quality threshold function call: function call
		Q1_lt_30 = Q30(quality_scores_R1)
		Q2_lt_30 = Q30(quality_scores_R2)
		qual_lt_30_249.extend((Q1_lt_30,Q2_lt_30))

		percent_reads_lt_30_R1 = Q1_lt_30/len(seqs1) * 100
		percent_reads_lt_30_R2 = Q2_lt_30/len(seqs2) * 100
		perc_qual_lt_30_249.extend((percent_reads_lt_30_R1,percent_reads_lt_30_R2))

		# Ambiguous base function call: function call
		countN1 = countN(seqs1)
		countN2 = countN(seqs2)
		ambi_calls_249.extend((countN1,countN2))
		
		# Descriptive stats for read1 length: function call
		r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read1_length)
		i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(read2_length)
		
		N_mean_249.extend((r_mean,i_mean))
		SD_249.extend((r_stdDev,i_stdDev))
		Variance_249.extend((r_var,i_var))
		median_249.extend((r_median,i_median))
		Q1_249.extend((r_Q1,i_Q1))
		Q3_249.extend((r_Q3,i_Q3))
		lwhisker_249.extend((r_lwhisker,i_lwhisker))
		hwhisker_249.extend((r_hwhisker,i_hwhisker))
		Skew_249.extend((r_skew,i_skew))
		G_mean_249.extend((r_gmean,i_gmean))

		
		# Descriptive stats for Q1 quality: function call
		q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores_R1)
		s_mean,s_stdDev,s_var,s_Q1,s_median,s_Q3,s_skew,s_gmean,s_lwhisker,s_hwhisker = stats(quality_scores_R2)
		
		qual_N_mean_249.extend((q_mean,s_mean))
		qual_SD_249.extend((q_stdDev,s_stdDev))
		qual_Variance_249.extend((q_var,s_var))
		qual_median_249.extend((q_median,s_median))
		qual_Q1_249.extend((q_Q1,s_Q1))
		qual_Q3_249.extend((q_Q3,s_Q3))
		qual_lwhisker_249.extend((q_lwhisker,s_lwhisker))
		qual_hwhisker_249.extend((q_hwhisker,s_hwhisker))
		qual_skew_249.extend((q_skew,s_skew))
		qual_G_mean_249.extend((q_gmean,s_gmean))
		
		
		# read Length thresholds: function call
		R1_le_149,R1_gt_149,R1_le_249,R1_gt_249,R1_le_299,R1_gt_299,tot_len1_le_149,tot_len1_gt_149,tot_len1_le_249,tot_len1_gt_249,tot_len1_le_299,tot_len1_gt_299 = threshold(read1_length)
		R2_le_149,R2_gt_149,R2_le_249,R2_gt_249,R2_le_299,R2_gt_299,tot_len2_le_149,tot_len2_gt_149,tot_len2_le_249,tot_len2_gt_249,tot_len2_le_299,tot_len2_gt_299 = threshold(read2_length)
	
		
		R_le_249.extend((R1_le_249,R2_le_249))
		R_gt_249.extend((R1_gt_249,R2_gt_249))
		
		
		
		
		perc_R1_le_249 = (R1_le_249/read_count1) * 100
		perc_R1_gt_249 = (R1_gt_249/read_count1) * 100
		perc_R2_le_249 = (R2_le_249/read_count2) * 100
		perc_R2_gt_249 = (R2_gt_249/read_count2) * 100
				
		final_perc_R1_le_249.extend((perc_R1_le_249,perc_R2_le_249))
		final_perc_R1_gt_249.extend((perc_R1_gt_249,perc_R2_gt_249))

		
		# avg_qual_score: function call
		avg_quality_1 = avg_qual_score(read1_length,quality_scores_R1)
		avg_quality_2 = avg_qual_score(read2_length,quality_scores_R2)
		
		final_avg_quality_249.extend((avg_quality_1,avg_quality_2))

		avg_length_1 = (tot_len1_le_249 + tot_len1_gt_249) / (R1_le_249 + R1_gt_249)
		avg_length_2 = (tot_len2_le_249 + tot_len2_gt_249) / (R2_le_249 + R2_gt_249)
		
		final_avg_length_249.extend((avg_length_1,avg_length_2))
	
		
		
				
	else:
		files_299.extend((f1,f2))
		
		
		# command line arguments
		fastq1 = f1
		fastq2 = f2

		# Parsing fastq: function call
		seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
		seqs2,quals2 = parseFastq(fastq2)
			
		# total number of reads
		read_count1 = len(seqs1)
		read_count2 = len(seqs2)
		
		total_no_reads_299.extend((read_count1,read_count2))
		
		# average quality scores for each read: function call
		read1_length,quality_scores_R1 = qual_score(quals1)
		read2_length,quality_scores_R2 = qual_score(quals2)
		
		
		# quality threshold function call: function call
		Q1_lt_30 = Q30(quality_scores_R1)
		Q2_lt_30 = Q30(quality_scores_R2)
		qual_lt_30_299.extend((Q1_lt_30,Q2_lt_30))

		percent_reads_lt_30_R1 = Q1_lt_30/len(seqs1) * 100
		percent_reads_lt_30_R2 = Q2_lt_30/len(seqs2) * 100
		perc_qual_lt_30_299.extend((percent_reads_lt_30_R1,percent_reads_lt_30_R2))

		# Ambiguous base function call: function call
		countN1 = countN(seqs1)
		countN2 = countN(seqs2)
		ambi_calls_299.extend((countN1,countN2))
		
		# Descriptive stats for read1 length: function call
		r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read1_length)
		i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(read2_length)
		
		N_mean_299.extend((r_mean,i_mean))
		SD_299.extend((r_stdDev,i_stdDev))
		Variance_299.extend((r_var,i_var))
		median_299.extend((r_median,i_median))
		Q1_299.extend((r_Q1,i_Q1))
		Q3_299.extend((r_Q3,i_Q3))
		lwhisker_299.extend((r_lwhisker,i_lwhisker))
		hwhisker_299.extend((r_hwhisker,i_hwhisker))
		Skew_299.extend((r_skew,i_skew))
		G_mean_299.extend((r_gmean,i_gmean))

		
		# Descriptive stats for Q1 quality: function call
		q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores_R1)
		s_mean,s_stdDev,s_var,s_Q1,s_median,s_Q3,s_skew,s_gmean,s_lwhisker,s_hwhisker = stats(quality_scores_R2)
		
		qual_N_mean_299.extend((q_mean,s_mean))
		qual_SD_299.extend((q_stdDev,s_stdDev))
		qual_Variance_299.extend((q_var,s_var))
		qual_median_299.extend((q_median,s_median))
		qual_Q1_299.extend((q_Q1,s_Q1))
		qual_Q3_299.extend((q_Q3,s_Q3))
		qual_lwhisker_299.extend((q_lwhisker,s_lwhisker))
		qual_hwhisker_299.extend((q_hwhisker,s_hwhisker))
		qual_skew_299.extend((q_skew,s_skew))
		qual_G_mean_299.extend((q_gmean,s_gmean))
		
		
		
		# read Length thresholds: function call
		R1_le_149,R1_gt_149,R1_le_249,R1_gt_249,R1_le_299,R1_gt_299,tot_len1_le_149,tot_len1_gt_149,tot_len1_le_249,tot_len1_gt_249,tot_len1_le_299,tot_len1_gt_299 = threshold(read1_length)
		R2_le_149,R2_gt_149,R2_le_249,R2_gt_249,R2_le_299,R2_gt_299,tot_len2_le_149,tot_len2_gt_149,tot_len2_le_249,tot_len2_gt_249,tot_len2_le_299,tot_len2_gt_299 = threshold(read2_length)
		
		
		R_le_299.extend((R1_le_299,R2_le_299))
		R_gt_299.extend((R1_gt_299,R2_gt_299))
		
		
		perc_R1_le_299 = (R1_le_299/read_count1) * 100
		perc_R1_gt_299 = (R1_gt_299/read_count1) * 100
		perc_R2_le_299 = (R2_le_299/read_count2) * 100
		perc_R2_gt_299 = (R2_gt_299/read_count2) * 100
		
		final_perc_R1_le_299.extend((perc_R1_le_299,perc_R2_le_299))
		final_perc_R1_gt_299.extend((perc_R1_gt_299,perc_R2_gt_299))
		
		#header.append("\n\n-----Stats for 299 bucket---------")
		
		# avg_qual_score: function call
		avg_quality_1 = avg_qual_score(read1_length,quality_scores_R1)
		avg_quality_2 = avg_qual_score(read2_length,quality_scores_R2)
		
		final_avg_quality_299.extend((avg_quality_1,avg_quality_2))

		avg_length_1 = (tot_len1_le_299 + tot_len1_gt_299) / (R1_le_299 + R1_gt_299)
		avg_length_2 = (tot_len2_le_299 + tot_len2_gt_299) / (R2_le_299 + R2_gt_299)
		
		final_avg_length_299.extend((avg_length_1,avg_length_2))

		
		
		
#function call
print_150bp()
print_250bp()
print_300bp()
