#!/usr/bin/env python3
# -*- coding: utf-8 -*-

####### SIDEROSCANNER #######
__author__ = 'Tom Stanton (tomdstanton@gmail.com)'
__version__ = '0.1'
__date__ = '31.03.20'

import os
import subprocess
from Bio import SeqIO
import argparse, textwrap
from argparse import RawTextHelpFormatter

### Args ###

def parse_args():
    
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description = textwrap.dedent('''\
        SideroScannner: A tool for accuratley annotating siderophore uptake proteins in bacteria.

        By Tom Stanton - Schneiders Lab - University of Edinburgh
        For queries/help/advice/comments/collaborations/chat:
        Email: T.D.Stanton@sms.ed.ac.uk --- Github: https://github.com/tomdstanton'''))
    
    parser.add_argument('--makedb', action = 'store_true', default = False,
                        help=('Setup TBDT DB and HMM profile'))
    
    parser.add_argument('-s', '--seq', nargs='*', type = str,#default = 'prot/SGH10.faa',
                        help=('Path to fasta input, autodetects DNA or Protein.'))
    
    parser.add_argument('-o', '--out', default='stdout', type = str,
                        help=('Path to output (comma-separated), otherwise prints to STDOUT.'))
    
    parser.add_argument('-d', '--db', default='irompdb', type = str,
                        help=('''Path to protein database:
If used with --makedb, will process your own protein fasta file into a compatible DB.
If used with scan, a diamond-formatted database is required.'''))
                        
    parser.add_argument('-m','--hmm', default='iromp.hmm', type = str,
                        help=('Path to HMM profile, must be pressed in hmmer format.'))
    
    parser.add_argument('-x', '--blastx', action = 'store_true', default = False,
                        help=('''Perform translated alignment (protein files will default to blastp).
WARNING: This is the fastest option for nucleotide files but will result in fewer hits.'''))
                              
    parser.add_argument('-b', '--blast_only', action = 'store_true', default = False,
                        help=('Turns off HMM pre-filter.'))
    
    parser.add_argument('-a','--annot', action = 'store_true', default = False,
                        help=('Add descriptions to hits'))
    
    parser.add_argument('-g', '--genloc', action = 'store_true', default = False,
                        help=('''Turns on chromosome/plasmid detection: 
Download the plsdb (https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)
Place the blast db files in sideroscanner script directory.
    -plsdb.fna.nsq
    -plsdb.fna.nhr
    -plsdb.fna.nin'''))				
    
    parser.add_argument('--logo', action = 'store_true', default = False)

    return parser.parse_args()
def main():
    
    makedb = parse_args().makedb
    db = parse_args().db
    out = parse_args().out
    seqs = parse_args().seq
    hmm = parse_args().hmm
    blast_only = parse_args().blast_only
    blastx = parse_args().blastx
    genloc = parse_args().genloc
    annot = parse_args().annot
    logo = parse_args().logo
    
    def print_logo():
        print('''
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████▌         ████████████████████████████████████████████          ████████
        ████████    ▄▄▄▄▄▄█████████████████████████████████████████████▄▄▄▄▄▄   ████████
        ████████   ▐█████████████████████████████▀███████████████████████████   ████████
        ████████   ▐████████▄█▀██████▐███▓███████▄██▀█████▌███▀██████████████   ████████
        █████████▄▄███████████▓▀█████▓██▌██████▌▀██▌▓█████▐███▓████████▀▌████▄▄▄████████
        ████████████████▌▀▀&▄Å▄B▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ÉÉÉÉBB▄Å$¢╬▀█████████████████
        ███████████████▄$ÉB▄▄██████████████████████████████████████▄▄É8$Ä█▓█████████████
        ██████████▀██▌╬Æ╬▄█████████████████████████████████████████████▄@Ñ╬█████████████
        ████████████▐Z╬▄████████████████████████████████████████████████▌ÑBZ█▀▀█████████
        ███████████▌ù⌂I██████████████████████████████████████████████████BÉ!w███████████
        ████                                                                        ████
        ████▄▄▄▄▄▄▄                                                          ╓▄▄▄▄▄▄████
        ████████████                                                         ███████████
        ███████████╛,                                                      ,▌ ██████████
        ███████████▄██▄                                                   '██▌▄█████████
        ████████████████ ,                                              ██ ▀████████████
        █████████▀▀████▌▄██ ▄▄▄▄ .,,,,, ¬▄▄▄▄▄ ,,,,  ▄▄▄  ▄▄▄  ▄▄▄▄ ▐██▌ █▄▐██▀█████████
        ████████   ▐█████,,▄█████▄█████▄▐█████ █████▄███▌ ███▌ ████▄ ████████   ████████
        ████████   ▐█████████████████████████████████████████████████████████   ████████
        ████████    ████████████████████████████████████████████████████████▀   ████████
        ████████          ████████████████████████████████████████████          ████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ███████  ~═~√^ ▀▀▀╙ █▀▀▀▀█▀▀▀▀█▀▀▀▀█  ══=¬█▀▀▀▀█▀▀▀▀█▀▀▀▀▀▀▀▀╙▀█▀▀▀▀▀▀▀▀▀███████
        ███████`"═²      ▀▀ ▄ ═∞═▐ ▐██▌ ▀▀  ``══  ▌ ▀▀"▀ =∞ ▐ ▐█▌ P ██ ▐ -∞══  █████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████
        ████████████████████████████████████████████████████████████████████████████████

        https://the-dots.com/users/deji-feyisetan-91484
        ''')

    FNULL = open(os.devnull, 'w') # This prevents stdout.

    if logo == True:
        if None not in (seqs, hmm, blast_only, blastx, genloc, out, annot, makedb):
            print('''--logo cannot be used with other arguments''')
        else:
            print_logo()
            exit()

    if makedb == True:
        if None not in (seqs, hmm, blast_only, blastx, genloc, out, annot):
            print('''--makedb cannot be used with other arguments except
            for creating your own DB with --db''')
        else:
            def make_iromp_db():
                def fetch(url, params, out):
                    r = requests.get(url,params)
                    open(out, 'wb').write(r.content)
                ### Create and process IROMP fasta file
                sequences = []   
                irompnames = ["iutA","fecA","fepA","cirA","iroN","fyuA","fhuA","fhuE","fiu","fatA","piuA","fptA",
                "pirA","fiuA","hxuC","bauA","pfeA","femA","foxA", "fitA","hmuR","pfuA","oprC","fpvA","fpvB"]
                for iromp in irompnames:
                    params = {'query': 'gene:{} AND taxonomy:Gammaproteobacteria AND ' \
                        'length:[500 TO 900] AND ' \
                        'locations:(location:"Cell outer membrane [SL-0040]")'.format(iromp),
                        'format': 'fasta'}
                    print("Fetching: "+iromp)
                    fetch("http://www.uniprot.org/uniprot/", params, iromp)
                    with open(iromp, "r") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            record.id = iromp
                            record.description = record.description.replace(record.name,'').lstrip()
                            sequences.append(record)
                SeqIO.write(sequences,'irompdb','fasta')
                for iromp in irompnames: 
                    os.remove(iromp)
                ### Get TBDT HMM from PFAM
                print("Fetching HMM")
                fetch('https://pfam.xfam.org/family/PF00593/hmm','','iromp.hmm')
                cmd = ['hmmpress','iromp.hmm']
                print("Pressing HMM")
                subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)

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

            if db == 'irompdb':
                import requests
                make_iromp_db()
                make_nr_db(db,db)
                os.remove(db)
            else:
                make_nr_db(db,db)
            
    if makedb == False:
        if seqs is None:
            print("Please provide at least one sequence")
            exit()
        import pandas as pd
        from io import StringIO
        
        def run_prodigal(in_file, out_file):
            cmd = ['prodigal','-i', in_file,'-a', out_file, '-q']
            print("Extracting proteins...")
            subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
            
        def run_hmmscan(in_file, hmm):
            from Bio import SearchIO
            cmd = ['hmmscan','--domtblout','hmmer_out', hmm, in_file]
            print("Filtering with HMM...")
            subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
            accessions = []
            for hit in SearchIO.parse('hmmer_out','hmmscan3-domtab'):
                accessions.append(hit._id)
            records = (r for r in SeqIO.parse(in_file,"fasta")if r.id in accessions)        
            SeqIO.write(records,'blast_in', "fasta")

        def run_diamond(method, in_file, db, filename, p):
            if annot == True:
                cmd = ['diamond', method, '-k', '1', '--id', '80', '--outfmt', '6', 'qseqid', 'salltitles',
                    'pident', 'bitscore', '--subject-cover', '40', '-q', in_file, '-d', db]

            if annot == False:
                cmd = ['diamond', method, '-k', '1', '--id', '80', '--outfmt', '6', 'qseqid', 'sseqid',
                    'pident', 'bitscore', '--subject-cover', '40', '-q', in_file, '-d', db]

            blast = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=FNULL)
            print("Aligning with "+method+"...")                
            blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
            diamond_df = pd.read_csv(blast_out, sep="\t", header=None,names=['Query','Hit','Percent_ID','Bitscore'])
            diamond_df.insert(0, 'Accession', filename)
            if genloc == True:
                if len(p) != 0: 
                    diamond_df['Genomic_Location'] = diamond_df['Query'].str.contains('|'.join(p))
                    diamond_df['Genomic_Location'] = diamond_df['Genomic_Location'].replace({True:'Plasmid',False:'Chromosome'})
                else:      
                    diamond_df['Genomic_Location'] = 'Chromosome'
            return(diamond_df)

        def run_screen(in_file):
            import shlex
            p = str(os.cpu_count())
            cmd = f'blastn -query {in_file} -task megablast -db plsdb.fna -outfmt \'6 qseqid\' ' \
            f'-num_threads {p} -perc_identity 90 -max_target_seqs 5 -qcov_hsp_perc 30'
            print("Screening for plasmids...")
            screen = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=FNULL)
            plasmids = screen.stdout.read().decode("utf-8").splitlines()
            plasmids = list(dict.fromkeys(plasmids))
            plasmid_number = len(plasmids)
            if plasmid_number == 0:
                print("No plasmids found")
            return plasmids

        ### Run the scan ###
        if out != 'stdout':
            df = pd.DataFrame()
        
        for seq in seqs:
            # Test to see if DNA or AA input
            test = next(SeqIO.parse(seq,"fasta"))
            filename = os.path.splitext(os.path.basename(seq))[0]

            # Protein input
            if "E" in test._seq._data:
                print("Protein fasta detected: "+filename)
                if genloc == True:
                    print("Cannot screen plasmids in protein fasta")
                else:
                    if blast_only == True:
                        x = run_diamond('blastp', seq, db, filename,0)
                        if out == 'stdout':
                            print(x.to_csv(sep='\t', index=False))
                        else:
                            df = df.append(x)
                    if blast_only == False:                   
                        run_hmmscan(seq, hmm)
                        x = run_diamond('blastp', 'blast_in', db, filename,0)
                        if out == 'stdout':
                            print(x.to_csv(sep='\t', index=False))
                        else:
                            df = df.append(x)
                
            # DNA input
            if "E" not in test._seq._data:
                print("Nucleotide fasta detected: "+filename)
                if blastx == True:
                    if genloc == True:
                        plasmids = run_screen(seq)
                    else:
                        plasmids = ""
                    x = run_diamond('blastx', seq, db, filename, plasmids)
                    if out == 'stdout':
                        print(x.to_csv(sep='\t', index=False))
                    else:
                        df = df.append(x)
                else:
                    seqname = os.path.splitext(seq)[0]
                    if blast_only == True:
                        if genloc == True:
                            plasmids = run_screen(seq)
                        else:
                            plasmids = ""
                        run_prodigal(seq, seqname+'.faa')
                        x = run_diamond('blastp', seqname+'.faa', db, filename, plasmids)
                        if out == 'stdout':
                            print(x.to_csv(sep='\t', index=False))
                        else:
                            df = df.append(x)
                    if blast_only == False:
                        if genloc == True:
                            plasmids = run_screen(seq)
                        else:
                            plasmids = ""
                        run_prodigal(seq, seqname+'.faa')
                        run_hmmscan(seqname+'.faa', hmm)
                        x = run_diamond('blastp', 'blast_in', db, filename, plasmids)
                        if out == 'stdout':
                            print(x.to_csv(sep='\t', index=False))
                        else:
                            df = df.append(x)
        if out != 'stdout':
            df.to_csv(out, index=False)
            print("Done, written to: "+out)
                
        ### Cleanup ###
        junk = ['hmmer_out','blast_in']
        for p in junk:
            if os.path.exists(p):
                os.remove(p)

if __name__ == "__main__":
    main()