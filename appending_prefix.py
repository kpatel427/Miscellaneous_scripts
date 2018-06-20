#!/usr/local/bin/python
# Khushbu Patel | 06/20/2018
#|__This script appends WGS_ids as a prefix to SRRs
#|__Load module Python/3.4
# create an excel file with two columns, FSIS IDs and corresponding SRRs and keep the column names as WGS_id and LabID, respectively.
# Usage: ./linelist.python

import pandas as pd
import subprocess

wgs_ID = []
SRR = []


data = pd.read_excel("FSIS.xlsx")

wgs_ID = data['WGS_id'].values.tolist()
SRR = data['LabID'].values.tolist()

for a,b in zip(wgs_ID,SRR):
	b = b+'.fastq.gz'	# Make sure you have the correct file extension
	subprocess.call(["mv",b,a+'_'+b])