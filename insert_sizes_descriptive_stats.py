#!/usr/local/bin/python
# Khushbu Patel | 05/29/2018
#|__This script requires samtools/1.4.1, Python 3.4 and installed python modules - numpy & scipy
#|__Converts sorted bam file to sam file format, extracts insert sizes of correctly paired mapped reads
#|__Provides descriptive stats for insert sizes 
# Usage: ./extract_quality_and_stats_fastq.py [input_sorted.bam]

import numpy as np
import sys
import subprocess
from scipy.stats import skew,mstats,uniform


inserts = []
insert_sizes = []

input_bam = sys.argv[1]
output_sam = "output.sam"
	
# To parse input BAM file
def parseBam(input_bam):
	subprocess.call(["samtools", "view", "-h", input_bam],stdout=open(output_sam,'w'))
	with open(output_sam,"r") as f:
		for line in f:
			if line.startswith('@'):	# skip the headers
				continue
			else:
				line = line.rstrip()
				temp = line.split()
				flag = temp[1]
				if(int(flag) == 99 or int(flag) == 163):
					insert_sizes.append(temp[8])
		return insert_sizes
	
# function call
inserts = parseBam(input_bam)


# To calculate descriptive stats
def stats(in_array):
	a = np.array(in_array)
	a = a.astype(int)	
	global min,max		# Gives error: local variable 'min' referenced before assignment
	min = min(a)
	max = max(a)
	
	mean = a.mean()
	std_dev = a.std()
	variance = np.var(a)
	Q1 = np.percentile(a,25)
	median = np.percentile(a,50)
	Q3 = np.percentile(a,75)
	skewness = skew(a)
	geometric_mean = mstats.gmean(a)
	uni_mean = uniform.mean(a)
	
	
	high = []
	low = []
	IQR = Q3 - Q1
	lower = Q1 - (1.5*IQR)
	upper = Q3 - (1.5*IQR)
	
	if(min < lower):
		low_whisker = min
	else:
		low_whisker = min
	
	if(max > upper):
		high_whisker = max
	else:
		high_whisker = max
	
	return mean,std_dev,variance,Q1,median,Q3,skewness,geometric_mean,low_whisker,high_whisker,uni_mean
	
	


# Descriptive stats for Insert sizes
i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker,i_uni_mean = stats(inserts)

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
print("Uniform mean = ", i_uni_mean)



