#!/usr/local/bin/python
# Khushbu Patel | 05/25/2018
# enter the pattern (p) and Alignment (t) to search from; returns the positions where exact matches occur
# Usage: ./exactmatch.py 

def exactMatch(p,t):
	occurences = []
	for i in range(len(t) - len(p)+1):
		match = True
		for j in range(len(p)):
			if(t[i+j] != p[j]):
				match = False
				break
		if match:
			occurences.append(i)
	return occurences
	
result = exactMatch('AT','ATTGCAGTCGATCGATTGCGATAGA')
print(result)
