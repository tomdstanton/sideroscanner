#!/usr/bin/env python3
# -*- coding: utf-8 -*-

####### SIDEROSCANNER #######
__author__ = 'Tom Stanton (tomdstanton@gmail.com)'
__version__ = '0.1'
__date__ = ''

import os
import subprocess
from Bio import SeqIO
import argparse, textwrap
from argparse import RawTextHelpFormatter

### Args ###
parser = argparse.ArgumentParser(description = textwrap.dedent('''\
        SideroScannner: A tool for accuratley annotating siderophore uptake proteins in bacteria.

        By Tom Stanton - Schneiders Lab - University of Edinburgh
        For queries/help/advice/comments/collaborations/chat:
        Email: T.D.Stanton@sms.ed.ac.uk --- Github: https://github.com/tomdstanton'''), formatter_class=RawTextHelpFormatter)
parser.add_argument('-b', '--blast_only', action = 'store_true',
                    help=('Use this option to skip HMM Pre-filter.'), default=False)
parser.add_argument('--hmm', default='iromp.hmm',
                    help=('Path to HMM profile, must be in hmm format.'), type = str)
parser.add_argument('--seq', default='SGH10.fna', 
                    help=('Path to fasta input, autodetects DNA or Protein.'), type = str)
parser.add_argument('--db', default='irompdb',
                    help=('Path to protein database, must be in fasta format.'), type = str)
parser.add_argument('--makedb', action = 'store_true',
                    help=('Setup DB and HMM profile') ,default=False)
parser.add_argument('--desc', action = 'store_true', 
                    help=('Add descriptions to hits'), default=False)

args = parser.parse_args()
hmm = args.hmm
seq = args.seq
db = args.db
makedb = args.makedb
blast_only = args.blast_only
desc = args.desc

### Global Variables ###

FNULL = open(os.devnull, 'w') # This prevents stdout.

if makedb == True:
    ### Make IROMP Database ###
    import requests
    def uniprot_search(iromp): # Sets up the UniProt REST API search
            params = {'query': 'gene:{} AND taxonomy:Gammaproteobacteria AND ' \
                      'length:[500 TO 900] AND ' \
                      'locations:(location:"Cell outer membrane [SL-0040]")'.format(iromp), # Makes sure the right type of proteins are being searched
                     'format': 'fasta'} # Downloads in fasta format
            response = requests.get("http://www.uniprot.org/uniprot/", params)
            f.write(response.content) # Writes output to file
    
    if db == 'irompdb':
        sequences = []   
        irompnames = ["iutA","fecA","fepA","cirA","iroN","fyuA","fhuA","fhuE","fiu","fatA",
                    "piuA","fptA","pirA","fiuA","hxuC","bauA","pfeA","femA","foxA", "fitA",
                    "hmuR","pfuA","oprC","fpvA","fpvB"]
        for iromp in irompnames: # Loops over IROMPs and writes to separate fasta file
            with open(iromp, 'wb') as f:
                print("Fetching: "+iromp)
                uniprot_search(iromp)
            with open(iromp, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    record.id = iromp
                    record.description = record.description.replace(record.name,'').lstrip()
                    sequences.append(record)
        irompdb = SeqIO.write(sequences,'irompdb','fasta')
        for iromp in irompnames: 
            os.remove(iromp)
    
    def make_nr_db(in_file, out_file):
        cmd = ['cd-hit','-i',in_file,'-o','nr','-c','1','-T','0']
        print('Removing redundancy')
        subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
        if out_file.endswith('.fasta'):
            out_file = os.path.splitext(out_file)[0]
        cmd = ['diamond','makedb','--in','nr','-d', out_file]
        print('Building DB: '+out_file)
        subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
        os.remove('nr.clstr')
        os.remove('nr')
    
    make_nr_db(db,db)
    if db == 'irompdb':
        os.remove(db)
    
    # Fetch HMM
    url = 'https://pfam.xfam.org/family/PF00593/hmm'
    print('Fetching HMM')
    r = requests.get(url, allow_redirects=True)
    open('iromp.hmm', 'wb').write(r.content)
    cmd = ['hmmpress','iromp.hmm']
    print('Processing HMM')
    press = subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
    print('Done')
        
if makedb == False:
    from io import StringIO
    from Bio import SearchIO
    import pandas as pd
    
    def run_prodigal(in_file, out_file):
        cmd = ['prodigal','-i', in_file,'-a', out_file, '-q']
        print(in_file+": DNA fasta detected, extracting proteins with Prodigal")
        subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
        
    def run_hmmscan(in_file, hmm):
        cmd = ['hmmscan','--domtblout','hmmer_out', hmm, in_file]
        print("Filtering with hmmscan")
        subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
        accessions = []
        for hit in SearchIO.parse('hmmer_out','hmmscan3-domtab'):
            accessions.append(hit._id)      
        records = (r for r in SeqIO.parse(in_file,"fasta")if r.id in accessions)        
        SeqIO.write(records,'blast_in', "fasta")

    def run_diamond(in_file, db):
        if desc == True:
            cmd = ['diamond', 'blastp', '-k', '1', '--id', '80', '--outfmt',
                   '6','qseqid','salltitles','pident','length','qstart','qend','sstart','send',
                   'evalue','bitscore','--subject-cover','80', '-q', in_file, '-d', db]
        if desc == False:
            cmd = ['diamond', 'blastp', '-k', '1', '--id', '80', '--outfmt',
                   '6','qseqid','sseqid','pident','length','qstart','qend','sstart','send',
                   'evalue','bitscore','--subject-cover','60', '-q', in_file, '-d', db]
        blast = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        print("Aligning proteins with Diamond")                  
        blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
        filename = os.path.splitext(seq)[0]
        df = pd.read_csv(blast_out, sep="\t", header=None, names=['Query','Hit','Percent ID','Length','Query Start','Query End','Subject Start','Subject End',
               'Evalue','Bitscore'])
        df.insert(0, 'Accession', filename)
        df.to_csv(filename+'_sideroscanner.csv', index = False)
        print("Done! File written to: "+filename+'_sideroscanner.csv')
        

    ### Run the scan ###
    # Test to see if DNA or AA input
    test = next(SeqIO.parse(seq,"fasta"))
    
    # Protein input with diamond only
    if "E" in test._seq._data and blast_only == True:
        print(seq+": Protein fasta detected")
        run_diamond(seq,db)
    
    # Protein input with hmmer and diamond
    if "E" in test._seq._data and blast_only == False:
        print(seq+": Protein fasta detected")
        run_hmmscan(seq, hmm)
        run_diamond('blast_in',db)
        
    # DNA input with diamond only
    if "E" not in test._seq._data and blast_only == True:
        seqname = os.path.splitext(seq)[0]
        run_prodigal(seq, seqname+'.faa')
        run_diamond(seqname+'.faa',db)
    
    # DNA input with hmmer and diamond
    if "E" not in test._seq._data and blast_only == False:
        seqname = os.path.splitext(seq)[0]
        run_prodigal(seq, seqname+'.faa')
        run_hmmscan(seqname+'.faa', hmm)
        run_diamond('blast_in',db)
    
    ### Cleanup ###
    junk = ['hmmer_out','blast_in']
    for p in junk:
        if os.path.exists(p):
            os.remove(p)
