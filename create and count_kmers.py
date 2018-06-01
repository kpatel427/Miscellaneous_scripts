#!/usr/local/bin/python
# Khushbu Patel | 06/1/2018
# Generates k-mer from a dataset of specific length, lexicographically sorted and provides count for occurences of a k-mer

import sys


def kmer_compo(dataset,num):

	result=[]
	count = {}
	if(num > len(dataset)):
		print("Second parameter cannot exceed the length of the dataset!")
		exit()
	for i in range(0,len(dataset)-num+1):
		result.append(dataset[i:i+num])
	result = sorted(result)
	
	for k in result:
		cnt = result.count(k)
		count[k] = cnt
		
	print(count)
	
kmer_compo(sys.argv[1],int(sys.argv[2]))