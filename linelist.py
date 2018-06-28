#!/usr/local/bin/python
#|__Khushbu Patel
#|__This script appends cluster code as a prefix to WGS_ids
#|__Load module Python/3.4
# create an excel file with two columns, FSIS IDs and corresponding SRRs and keep the column names as WGS_id and LabID, respectively.
# Usage: ./linelist.python

import pandas as pd
import subprocess
import glob
import os

temp = []
wgs_ID = []
SRR = []
files = []
basename = []

for f in glob.glob('*analreq*.xlsx'):
	files.append(f)
	data = pd.read_excel(f)
	temp = data['WGS ID'].values.tolist()
	basename.append(os.path.splitext(f)[0])

for a in basename:
	for x in temp:
		if x.startswith('PNU'):
			print(a,x)
		elif x.startswith('FSIS'):
			print(a,x)
		else:
			continue
