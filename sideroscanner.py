#!/usr/bin/env python3
# -*- coding: utf-8 -*-

####### SIDEROSCANNER #######
__author__ = 'Tom Stanton (tomdstanton@gmail.com)'
__version__ = '0.1'
__date__ = '15.04.20'

import os, shlex, requests, argparse, textwrap, subprocess, re, sys
from Bio import SeqIO, SearchIO, Entrez
from argparse import RawTextHelpFormatter
import pandas as pd
from io import StringIO
#from Bio.Align.Applications import MuscleCommandline as muscle

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description = textwrap.dedent('''\
        SideroScannner: A tool for accuratley annotating siderophore uptake proteins in bacteria.

        By Tom Stanton - Schneiders Lab - University of Edinburgh
        For queries/feedback/comments/collaborations/chat/coffee donation:
        Email: T.D.Stanton@sms.ed.ac.uk --- Github: https://github.com/tomdstanton
        
        Please note:
            Optimising the semi-curated protein database is an ongoing process.'''))
    
    parser.add_argument('--makedb', action = 'store_true', default = False,
                        help=('Setup TBDT DB and HMM profile'))
    
    parser.add_argument('-s', '--seq', nargs='*', type = str,
                        help=('Path to fasta input, detects DNA or Protein.'))
    
    parser.add_argument('-o', '--out', default='stdout', type = str,
                        help=('Path to output (comma-separated), otherwise prints to STDOUT.'))
    
    parser.add_argument('-d', '--db', default='irompdb', type = str,
                        help=('''Path to protein database:
If used with --makedb, will process your own protein fasta file into a compatible DB.
If used with scan, a diamond-formatted database is required.'''))
                        
    parser.add_argument('-m','--hmm', default='iromp.hmm', type = str,
                        help=('Path to HMM profile, must be "pressed" in hmmer format.'))
    
    parser.add_argument('-x', '--blastx', action = 'store_true', default = False,
                        help=('''Perform translated alignment (protein files will default back to blastp).
WARNING: This is the fastest option for nucleotide files but will result in fewer hits.'''))
                              
    parser.add_argument('-b', '--blast_only', action = 'store_true', default = False,
                        help=('Turns off HMM pre-filter.'))
    
    parser.add_argument('-a','--annot', action = 'store_true', default = False,
                        help=('Add descriptions to hits.'))
    
    parser.add_argument('-f','--fur', action = 'store_true', default = False,
                        help=('''If using contigs/assemblies, you can scan the
promoter regions of hits for Fur binding sites.'''))
    
    parser.add_argument('-p','--pwm', default='pwm', type = str,
                        help=('''Path to meme formatted position weight matrix.
Use with --fur'''))
    
    parser.add_argument('-g', '--genloc', action = 'store_true', default = False,
                        help=('''Turns on chromosome/plasmid detection: 
Download the plsdb (https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)
Place the blast db files in sideroscanner script directory.
    -plsdb.fna.nsq
    -plsdb.fna.nhr
    -plsdb.fna.nin'''))				
    
    parser.add_argument('--logo', action = 'store_true', default = False)
    if len(sys.argv)==1:
        parser.print_help()
        print('\n'+'^^^ No flags given, see help above ^^^')
        sys.exit(1)
    return parser.parse_args()

### Globals ###
makedb = parse_args().makedb
db = parse_args().db
out = parse_args().out
seqs = parse_args().seq
hmm = parse_args().hmm
blast_only = parse_args().blast_only
blastx = parse_args().blastx
genloc = parse_args().genloc
annot = parse_args().annot
fur = parse_args().fur
pwm = parse_args().pwm
logo = parse_args().logo
FNULL = open(os.devnull, 'w') # This prevents stdout.
p = str(os.cpu_count())
   
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
           
def fetch(url, params):
    q = requests.get(url,params)
    return q.content

def make_iromp_db():
    sequences = []
    irompnames = ["iutA","fecA","fepA","cirA","iroN","fyuA","fhuA","fhuE","fiu","fatA","piuA","fptA",
    "pirA","fiuA","hxuC","bauA","pfeA","femA","foxA", "fitA","hmuR","pfuA","oprC","fpvA","fpvB",
    "fcuA","btuB"]
    for iromp in irompnames:
        
    ### Uniprot Proteins ###
        # params = {'query': 'gene_exact:{} AND taxonomy:Gammaproteobacteria AND ' \
        #     'length:[500 TO 900] AND locations:(location:"Cell outer membrane [SL-0040]")'.format(iromp,iromp),
        #     'format': 'fasta'}
        # print("Fetching: "+iromp)
        # handle1 = fetch("http://www.uniprot.org/uniprot/", params).decode('utf-8').splitlines()
        # records = list(SeqIO.parse(handle1, "fasta"))
        # for r in records:
        #     r.id = iromp
        #     r.description = r.description.replace(r.name,'').lstrip()
        #     sequences.append(r)

    ### Entrez Proteins ###
        Entrez.email = "tomdstanton@gmail.com"
        handle=Entrez.esearch(db="ipg",
        term = '{}[Protein Name] AND prokaryotes[filter]' \
        '500:900[Sequence Length]'.format(iromp),
        retmax = 10000000, idtype="acc", usehistory="y")
        protein_id=Entrez.read(handle)['IdList']
        for prid in protein_id:
            result = Entrez.efetch(db='ipg', id=prid, retmax = 10000000, rettype='fasta', retmode='text')
            record = SeqIO.read(result, 'fasta')
            record.id = iromp
            sequences.append(record)
    return SeqIO.write(sequences,'irompdb','fasta') 
                
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
                    
def run_prodigal(in_file, out_file):
    cmd = ['prodigal','-i', in_file,'-a', out_file, '-q']
    print("Extracting proteins...")
    return subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)

def run_mast(hits, seq, in_file):
    print("Extracting blast hit promoter sequences...")
    blast_hits = hits['Query'].tolist()
    records = list(SeqIO.parse(in_file,"fasta"))
    bed_file = pd.DataFrame(columns = ["chrom", "start",
                                       "end", "strand"])
    for r in records:    
        if r.id in blast_hits:
            headers = r.description + "\n"
            bed = re.split('\#|\s',headers.replace(' ',''))[0:4]
            a_series = pd.Series(bed, index = bed_file.columns)
            bed_file = bed_file.append(a_series, ignore_index=True)
    bed_file['chrom'] = bed_file['chrom']
    bed_file.loc[bed_file['strand'] == '1',
                 'flank_start'] = bed_file['start'].astype(int) - 450
    bed_file.loc[bed_file['strand'] == '1',
                 'flank_end'] = bed_file['start'].astype(int) + 150
    bed_file.loc[bed_file['strand'] == '-1',
                 'flank_start'] = bed_file['end'].astype(int) -150
    bed_file.loc[bed_file['strand'] == '-1',
                 'flank_end'] = bed_file['end'].astype(int) + 450
    chrs = {}
    for r in SeqIO.parse(seq, "fasta"):
        chrs[r.id] = r.seq
    promoter_reg = ''
    for index, row in bed_file.iterrows():
        hit_accession = row['chrom']
        chracc = hit_accession.split('_')
        b = "_".join(chracc[:2])
        f_start = int(row['flank_start'])
        f_end = int(row['flank_end'])
        flank = chrs[b][f_start:f_end]
        promoter_reg = promoter_reg+ \
            '\n'+'>'+hit_accession+':'+str(f_start) \
                +'-'+str(f_end)+'\n'+str(flank)
    mast_in = promoter_reg.format("fasta").strip()
    print('Performing motif alignment...')
    cmd = ['mast', '-hit_list', pwm, '-']
    mast = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    mast.stdin.write(mast_in.encode())
    mast_out = StringIO(mast.communicate()[0].decode('utf-8'))   
    mast_df = pd.read_csv(mast_out, sep=" ", header=1)
    mast_df = mast_df[:-1]
    mast_df[['Query','flank']] = mast_df['#'].str.split(":",expand=True)
    mast_df[['flank_start','flank_end']] = mast_df['flank'].str.split("-",expand=True)
    mast_df = mast_df.drop(columns=['#','id','score', 'flank', 'hit_end'])
    mast_df = mast_df.rename(columns={"sequence_name":"strand+/-",
                            "(strand+/-)motif": "motif",
                            "alt_id":"motif_start",
                            "hit_start":"motif_end",
                            "hit_p-value":"motif_p-value"})
    return mast_df[['Query','flank_start','flank_end','motif',
                    'strand+/-','motif_start', 'motif_end',
                    'motif_p-value']]

def run_hmmsearch(in_file):
    cmd = ['hmmsearch', hmm, in_file]
    print("Filtering with HMM...")
    hmmsearch = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=FNULL)
    hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
    qresult = next(SearchIO.parse(hmmer_out,'hmmer3-text'))
    return qresult.hit_keys

def run_diamond(method, in_file, filename, accessions, plasmids):
    if annot == True:
        hit = 'salltitles'
    else:
        hit = 'sseqid'
    if blast_only == False:
        sequences = ''
        records = list(SeqIO.parse(in_file,"fasta"))
        for r in records:    
            if r.id in accessions:
                sequences = sequences + r.format("fasta")
        cmd = ['diamond', method, '-k', '1', '--id', '95', '--outfmt', '6', 'qseqid',
               hit, 'pident', 'bitscore', '--subject-cover', '40', '-d', db]
        blast = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,stderr=FNULL)
        blast.stdin.write(sequences.encode())
        blast_out = StringIO(blast.communicate()[0].decode('utf-8'))   
    else: 
         cmd = ['diamond', method, '-k', '1', '--id', '95', '--outfmt', '6', 'qseqid', hit,
                'pident', 'bitscore', '--subject-cover', '40', '-q', in_file, '-d', db]
         blast = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=FNULL)
         blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    print("Aligning with "+method+"...")                
    diamond_df = pd.read_csv(blast_out, sep="\t", header=None,names=['Query','Hit','Percent_ID','Bitscore'])
    diamond_df.insert(0, 'Accession', filename)
    if genloc == True:
        if len(plasmids) != 0: 
            diamond_df['Genomic_Location'] = diamond_df['Query'].str.contains('|'.join(plasmids))
            diamond_df['Genomic_Location'] = diamond_df['Genomic_Location'].replace({True:'Plasmid',False:'Chromosome'})
        else:      
            diamond_df['Genomic_Location'] = 'Chromosome'
    return diamond_df

def run_screen(in_file):
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

#def run_muscle():
            
def main():
    if logo == True:
        if None not in (seqs, hmm, blast_only, blastx, genloc, out, annot, makedb, db):
            print('''--logo cannot be used with other arguments''')
        else:
            print_logo()
            sys.exit(1)

    if makedb == True:
        if None not in (seqs, hmm, blast_only, blastx, genloc, out, annot, logo):
            print('''--makedb cannot be used with other arguments except
            for creating your own DB with --db''')
            sys.exit(1)
        if db == 'irompdb':
            make_iromp_db()
            make_nr_db(db,db)
            #os.remove(db)
        else:
            make_nr_db(db,db)
        print("Fetching HMM")
        open('iromp.hmm', 'wb').write(fetch('https://pfam.xfam.org/family/PF00593/hmm',''))
        cmd = ['hmmpress','iromp.hmm']
        print("Pressing HMM")
        subprocess.run(cmd, stdout=FNULL, stderr=subprocess.STDOUT)

    if makedb == False:
        if len(seqs) == 0:
            print("Please provide at least one sequence with -s")
            sys.exit(1)
        print('Using '+p+' cores...' )
        if out != 'stdout':
            df = pd.DataFrame()

        for seq in seqs:
            # Test to see if DNA or AA input
            test = next(SeqIO.parse(seq,"fasta"))
            filename = os.path.splitext(os.path.basename(seq))[0]
            accessions = ''
            plasmids = ''
            
            # Protein input
            if "E" in test._seq._data:
                print("Protein fasta detected: "+filename)                   
                if genloc == True:
                    print("Cannot screen plasmids in protein fasta")
                    sys.exit(1)
                if fur == True:
                    print("Cannot screen TFBS in protein fasta")
                    sys.exit(1) 
                if blast_only == False:                   
                    accessions = run_hmmsearch(seq)
                x = run_diamond('blastp', seq, filename, accessions, 0)
                
            # DNA input
            if "E" not in test._seq._data:
                print("Nucleotide fasta detected: "+filename)
                seqname = os.path.splitext(seq)[0]
                if genloc == True:
                    plasmids = run_screen(seq)
                if blastx == True:
                    x = run_diamond('blastx', seq, filename, accessions, plasmids)
                else:
                    run_prodigal(seq, seqname+'.faa')
                    if blast_only == False:
                        accessions = run_hmmsearch(seqname+'.faa')
                    x = run_diamond('blastp', seqname+'.faa', filename, accessions, plasmids)
                if fur == True:                    
                    x = pd.merge(x, run_mast(x, seq, seqname+'.faa'), on='Query')
 
            if out == 'stdout':
                print(x.to_csv(sep='\t', index=False))
            else:
                df = df.append(x)
        if out != 'stdout':
            df.to_csv(out, index=False)
            print("Done, written to: "+out)

if __name__ == "__main__":
    main()