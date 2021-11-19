#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#how to use: masking.py target.fasta blasthit.tsv
#blasthit.tsv is the file created by blast against database built from target.fasta with '-outfmt 6' option.

import os
import tempfile
import sys
import csv
import shutil
import pandas as pd
from Bio import SeqIO

fasta_in = sys.argv[1]
query = sys.argv[2]

masked = "./masked.fasta"
shutil.copy(fasta_in, masked)

df = pd.read_csv(query, usecols=[1, 8, 9], delimiter='\t', header=None)
df.to_csv('subject.tsv')

f = open('subject.tsv','r')
reader = csv.reader(f)
for row in reader:
	temp = tempfile.TemporaryFile(mode='w+t')
	sub = row[1]
	for record in SeqIO.parse(masked, 'fasta'):
		id_part = record.id
		m_part = id_part.rstrip()
		description_part = record.description
		seq = record.seq
		if m_part == sub:
			mutable_seq = seq.tomutable()
			if int(row[2]) < int(row[3]):
				start_raw = row[2]
				start = int(start_raw) - 1
				end = row[3]
			else:
				start_raw = row[3]
				start = int(start_raw) - 1
				end = row[2]
			for num in range(int(start), int(end)): #律速。replaceで速くなる。
				mutable_seq[num] = "N"
			fasta_seq = '>' + description_part + '\n' + mutable_seq + '\n'
		else:
			fasta_seq = '>' + description_part + '\n' + seq + '\n'
		temp.writelines(fasta_seq)
	temp.seek(0)
	masked = open("masked.fasta", "w")
	contents = temp.readlines()
	masked.writelines(contents)
	masked.close()
	masked = open("masked.fasta")
	temp.close()
f.close()