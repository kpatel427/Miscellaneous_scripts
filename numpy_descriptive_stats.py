#!/usr/local/bin/python
# Khushbu Patel | 05/29/2018
# numpy descriptive stats

import numpy as np

arr = [5,10,15,20,25]

def stats(in_array):
	a = np.array(in_array)
	mean = a.mean()
	std_dev = a.std()
	variance = np.var(a)
	Q1 = np.percentile(a,25)
	median = np.percentile(a,50)
	Q3 = np.percentile(a,75)
	
	return mean,std_dev,variance,Q1,median,Q3
	
	
results = stats(arr)
print(results)