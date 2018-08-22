#!/usr/bin/Python
# DE Analysis for Target MYCN data
# Finds common and unique column names from count-matrix and Target-MYCN-status file anmd creates a new file with unique col names with their MYCN status

lst = [] #target_mycn
status = []
lst2 = []
result = []
cnt = 0

with open("path/to/target_MYCN_status.txt",'r') as target:
	with open("path/to/unique_column_names.txt",'r') as colnames:
		for line in target:	
				col_1 = line.split("\t")
				lst.append(col_1[0])
				status.append(col_1[1])
		for col in colnames:
			col.rstrip("\n")
			lst2.append(col)


for a,x in zip(lst,status):
	for b in lst2:
		if(a in b):
			#print(a+" "+x)
			fh = open("path/to/New_target_MYCN_status.txt",'a')
			fh.write(a+" "+x)
