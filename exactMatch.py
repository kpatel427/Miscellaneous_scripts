#!/usr/local/bin/python
# Khushbu Patel | 05/25/2018
# enter the pattern (p) and Alignment (t) to search from; returns the positions where exact matches occur
# Usage: ./exactmatch.py 


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t



def exactMatch(p,t):
	occurences = []
	reverse = reverseComplement(p)
	for i in range(len(t) - len(p)+1):
		match = True
		for j in range(len(p)):
			if(t[i+j] != p[j]):
				match = False
				if(t[i+j] != reverse[j]):
					match = False
					break
				else:
					match = True
					
		if match:
			occurences.append(i)
	return occurences
