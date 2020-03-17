#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'Tom Stanton (tomdstanton@gmail.com)'
__version__ = '0.1'
__date__ = ''

import os
from Bio import SearchIO, SeqIO
from io import StringIO
import argparse
import subprocess
import pandas as pd

parser = argparse.ArgumentParser(description='SideroScanner: a tool for annotating siderophore uptake proteins in bacteria')
parser.add_argument('-p', '--threads', default='4', type = str)
parser.add_argument('--hmm', default='iromp.hmm', type = str)
parser.add_argument('--seq', default='TS-SGH10.faa', type = str)
parser.add_argument('--db', default='irompdb', type = str)
args = parser.parse_args()
threads = args.threads
hmm = args.hmm
seq = args.seq
irompdb = args.db

FNULL = open(os.devnull, 'w')
hmmscan_cmd = ['hmmscan','--cpu', threads,
               '--domtblout','hmmerout', hmm, seq]
run_hmmscan = subprocess.run(hmmscan_cmd, stdout=FNULL, stderr=subprocess.STDOUT)
accessions = []
for hit in SearchIO.parse('hmmerout', 'hmmscan3-domtab'):
    accessions.append(hit._id)      
records = (r for r in SeqIO.parse(seq,"fasta")if r.id in accessions)        
blastp_in = SeqIO.write(records, 'blast_in', "fasta")
blast_cmd = ['diamond','blastp','-p',
             threads,'-q', 'blast_in', '-d', irompdb]
run_blast = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE)
blast_out = StringIO(run_blast.communicate()[0].decode('utf-8'))
df = pd.read_csv(blast_out, sep="\t", header=None)
df_filtered = df[(df[2] >= 80.0)]
idx = df_filtered.groupby([0])[11].transform(max) == df_filtered[11]
out = df_filtered[idx].drop_duplicates([0])
print(out.to_csv(sep='\t', index=False,header=False))
os.remove('hmmerout')
os.remove('blast_in')
