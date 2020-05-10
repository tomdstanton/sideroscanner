#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton (tomdstanton@gmail.com)'
__version__ = '0.1'
__date__ = '01.05.20'

import os, requests, argparse, textwrap, subprocess, re, sys, warnings
from Bio import SeqIO, SearchIO, BiopythonWarning
from argparse import RawTextHelpFormatter
import pandas as pd
from io import StringIO
from IPython.display import display
from Bio.SeqUtils.ProtParam import ProteinAnalysis

warnings.simplefilter('ignore', BiopythonWarning)

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description = textwrap.dedent('''\
        SideroScannner: A tool for annotating IROMPs in bacteria.

        By Tom Stanton - Schneiders Lab - University of Edinburgh
        For queries/feedback/comments/collaborations/chat/coffee donation:
        Email: T.D.Stanton@sms.ed.ac.uk --- Github: https://github.com/tomdstanton
        
        Please note:
            Curating the IROMP HMM library is an ongoing process.'''))
    
    parser.add_argument('-t', default=0, type = int,
                        help=('Number of threads to use, will default to all'))
    
    parser.add_argument('-s', nargs='*', type = str,
                        help=('Fasta input, accepts multiple DNA or AA files.'))
    
    parser.add_argument('-o', default='stdout', type = str,
                        help=('Output file (comma-separated)'))
    
    parser.add_argument('-g', action = 'store_true', default = False,
                        help=('Turns on chromosome/plasmid detection'))
    
    parser.add_argument('-f', action = 'store_true', default = False,
                        help=('If using contigs/assemblies, scan the promoter regions of hits for Fur binding sites'))
    
    parser.add_argument('-w', default='', type = str,
                        help=('Write proteins from analysis to file. Use "seq" for default filename'))
    
    parser.add_argument('-l', action = 'store_true', default = False,
                        help=('Adjusts gene prediction and only plug HMM filter for low-quality assemblies'))
    
    parser.add_argument('-p','--pwm', default='pwm', type = str,
                        help=('Path to meme formatted position weight matrix (use with --fur)'))
    
    parser.add_argument('-i', default='iromps.hmm', type = str,
                        help=('IROMP HMM, must be "pressed" in hmmer3 format'))
                        
    parser.add_argument('-d', default=['PF07715.hmm', 'PF00593.hmm'],
                        type = list,  help=('hmmer3 formatted Plug and TonB dependent receptor domain HMMs'))			
    
    parser.add_argument('--logo', action = 'store_true', default = False)
    
    if len(sys.argv)==1:
        parser.print_help()
        print('\n'+'Please provide at least one argument or -h for help'+'\n')
        sys.exit(1)
    return parser.parse_args()

### Globals
threads = parse_args().t
out = parse_args().o
seqs = parse_args().s
iromps = parse_args().i
domains = parse_args().d
lowq = parse_args().l
genloc = parse_args().g
fur = parse_args().f
pwm = parse_args().pwm
logo = parse_args().logo
write = parse_args().w
   
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

def run_prodigal(infile, quality):
    print("Extracting proteins...")
    cmd = ['prodigal','-i', infile, '-o', '/dev/null', '-a',
           '/dev/stdout', '-q', '-p', quality]
    prodigal = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return prodigal.communicate()[0].decode('utf-8')

def run_hmmsearch(molecule, infile, cpus, quality):
    if molecule == 'aa':
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus, domains[0], infile]
        hmmsearch = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        qresult = next(SearchIO.parse(hmmer_out,'hmmer3-text'))
        records = list(SeqIO.parse(infile,"fasta"))
    else:
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus, domains[0], '-']
        hmmsearch = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
        hmmsearch.stdin.write(infile.encode())
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        qresult = next(SearchIO.parse(hmmer_out,'hmmer3-text'))
        records = SeqIO.parse(StringIO(infile), 'fasta')
    filter1 = ''
    for r in records:    
        if r.id in qresult.hit_keys:
            filter1 = filter1 + r.format("fasta")
    print("Filtered "+str(filter1.count('>'))+" proteins with "+domains[0])
    if quality == 'meta' or len(filter1) == 0:
        return filter1
    else:
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus, domains[1], '-']
        hmmsearch = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE)
        hmmsearch.stdin.write(filter1.encode())
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        qresult = next(SearchIO.parse(hmmer_out,'hmmer3-text'))
        records = SeqIO.parse(StringIO(filter1), 'fasta')        
        filter2 = ''
        for r in records:    
            if r.id in qresult.hit_keys:
                filter2 = filter2 + r.format("fasta")
        print("Filtered "+str(filter2.count('>'))+" proteins with "+domains[1])
        return filter2

def run_hmmscan(infile, cpus, molecule):
    print("Scanning against HMM library...")
    hmmscan_df = pd.DataFrame(columns = ['contig','query', 'hit', 'hit_range',
                                         'score','evalue'])
    cmd = ['hmmscan', '-E', '1e-50', '--cpu', cpus, iromps, '-']
    hmmscan = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)
    hmmscan.stdin.write(infile.encode())
    hmmer_out = StringIO(hmmscan.communicate()[0].decode('utf-8'))
    for qresult in SearchIO.parse(hmmer_out, "hmmer3-text"):
        if len(qresult.hits) > 0:
            hmmscan_df = hmmscan_df.append(
                {'contig':qresult._id.rsplit('_', 1)[0], 'query':qresult._id,
                 'hit':qresult.hits[0]._id, 'score':qresult.hits[0].bitscore,
                 'evalue':"%.3g" % qresult.hits[0].evalue, 
                 'hit_range':str(qresult.hits[0].hsps[0].hit_range[0]) \
                     +'-'+str(qresult.hits[0].hsps[0].hit_range[1])},
                ignore_index = True)
    if hmmscan_df.empty == True:
        print("No proteins annotated")
        return hmmscan_df
    
    weight_df = pd.DataFrame(columns = ['query','mass(Da)'])
    for r in SeqIO.parse(StringIO(infile), 'fasta'):
        X = ProteinAnalysis((r._seq._data.replace('*', '')))
        weight_df = weight_df.append({
            'query':r.id, 'mass(Da)':int(X.molecular_weight())},
            ignore_index = True)
    hmmscan_df = hmmscan_df.merge(weight_df, on = 'query', how = 'left')
    if molecule == 'aa':
        return hmmscan_df
    else:    
        bed_df = pd.DataFrame(columns = ['query', 'start', 'end', 'strand'])
        for r in SeqIO.parse(StringIO(infile), 'fasta'):
            bed = re.split('\#|\s', (r.description + "\n").replace(' ',''))[0:4]
            bed_df = bed_df.append(pd.Series(
                bed,index=['query', 'start', 'end', 'strand']),
                ignore_index=True)
        return pd.merge(hmmscan_df, bed_df, on = 'query')

def plasmid_screen(infile, cpus, hits):
    print("Screening for plasmids...")
    cmd = ['blastn',
           '-query', infile,
           '-db', 'plsdb.fna',
           '-culling_limit', '1',
           '-outfmt', '6',
           '-max_hsps', '1',
           '-evalue', '1e-50',
           '-num_threads',  cpus, 
           '-perc_identity','90',
           '-qcov_hsp_perc', '10',
           '-max_target_seqs', '1']
    blast = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    plasmid_df = pd.DataFrame(columns = ['contig', 'plasmid/mge', 'span'])
    for qresult in SearchIO.parse(blast_out, "blast-tab"):
       for h in qresult.hits:
           plasmid_df = plasmid_df.append(
               {'contig':h.query_id,'plasmid/mge':h.id,
                'span':str(h.hsps[0].query_range[0])+'-'+str(h.hsps[0].query_range[1])},
               ignore_index = True)
    hits = hits.merge(plasmid_df, how='left')
    if plasmid_df.empty == True:
        print("No plasmids found")
    return hits.fillna(value={'plasmid/mge': 'Chromosome'})

def mge_screen(infile, cpus, plasmids):
    print("Screening for MGEs...")
    cmd = ['blastn', 
           '-query', infile,
           '-db', 'mge',
           '-culling_limit', '1', 
           '-outfmt', '6',
           '-max_hsps', '1',
           '-evalue', '1e-50',   
           '-num_threads', cpus,
           '-perc_identity', '90']
    blast = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    mge_df = pd.DataFrame(columns = ['contig', 'plasmid/mge', 'span'])
    for qresult in SearchIO.parse(blast_out, "blast-tab"):
       for h in qresult.hits:
           name = h.id.split('|')[2]
           mge_df = mge_df.append(
               {'contig':h.query_id,'plasmid/mge':name,
                'span':h.hsps[0].query_range}, ignore_index = True)
    chromosome_contigs = list(plasmids.loc[plasmids['plasmid/mge'] == 'Chromosome','contig'])
    mge_df = mge_df[mge_df['contig'].isin(chromosome_contigs)]
    if mge_df.empty == True:
        print("No MGEs found")
    else:
        def range_subset(range1, range2):
          if not range1:
              return True 
          if not range2:
              return False
          if len(range1) > 1 and range1.step % range2.step:
              return False
          return range1.start in range2 and range1[-1] in range2
        mge_dict = dict(zip(mge_df['plasmid/mge'], mge_df['span']))
        for row in plasmids.itertuples():
            for key in mge_dict:
                if range_subset(range(int(row.start), int(row.end)),
                                 range(mge_dict[key][0], mge_dict[key][1]))==True:
                    plasmids.at[row.Index, 'plasmid/mge'] = key
                    plasmids.at[row.Index, 'span'] = \
                        str(mge_dict[key][0])+'-'+str(mge_dict[key][1])
    return plasmids
            
def run_mast(hits, infile):
    print('Searching for motif hits in promoter regions...')
    hits.loc[hits['strand'] == '1',
                 'flank_start'] = hits['start'].astype(int) - 450
    hits.loc[hits['strand'] == '1',
                 'flank_end'] = hits['start'].astype(int) + 150
    hits.loc[hits['strand'] == '-1',
                 'flank_start'] = hits['end'].astype(int) -150
    hits.loc[hits['strand'] == '-1',
                 'flank_end'] = hits['end'].astype(int) + 450
    chrs = {}
    for r in SeqIO.parse(infile, "fasta"):
        chrs[r.id] = r.seq
    promoter_reg = ''
    for index, row in hits.iterrows():
        query = row['query']
        chracc = row['contig']
        f_start = int(row['flank_start'])
        f_end = int(row['flank_end'])
        flank = chrs[chracc][f_start:f_end]
        promoter_reg = promoter_reg+ \
            '\n'+'>'+query+':'+str(f_start) \
                +'-'+str(f_end)+'\n'+str(flank)
    promoters = promoter_reg.format("fasta").strip()
    cmd = ['mast', '-hit_list', pwm, '-']
    mast = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    mast.stdin.write(promoters.encode())
    mast_out = StringIO(mast.communicate()[0].decode('utf-8'))   
    mast_df = pd.read_csv(mast_out, sep=" ", header=1)
    if mast_df.empty == True:
        print("No TFBS found")
        return hits.drop(columns=['flank_start','flank_end'])
    else:
        mast_df = mast_df[:-1]
        mast_df[['query','flank']] = mast_df['#'].str.split(":",expand=True)
        mast_df[['flank_start','flank_end']] = mast_df['flank'].str.split("-",expand=True)
        mast_df = mast_df.drop(columns=['#','id','score', 'flank', 'hit_end', 'sequence_name'])
        mast_df = mast_df.rename(columns={"(strand+/-)motif": "motif",
                                "alt_id":"motif_start",
                                "hit_start":"motif_end",
                                "hit_p-value":"motif_p-value"})
        mast_df['motif_start'] = (mast_df['motif_start'].astype(int) + mast_df['flank_start'].astype(int)).astype(str)
        mast_df['motif_end'] = (mast_df['motif_end'].astype(int) + mast_df['flank_start'].astype(int)).astype(str)
        mast_df = mast_df.drop(columns=['flank_start','flank_end'])
        hits = hits.drop(columns=['flank_start','flank_end'])
        hits = hits.merge(mast_df, how = 'left', on = 'query')
        return hits

def grab_proteins(infile, hits, path, seq):
    to_write = []
    try:
        records = SeqIO.parse(StringIO(infile), 'fasta')
        for r in records:    
            write_df = hits[hits['query'].str.match(r.id)]
            r.id = write_df['query'].values[0]+'_'+write_df['hit'].values[0]
            r.description = ''
            to_write.append(r)
        if path == "seq":
            path = str(os.path.splitext(os.path.basename(seq))[0])+"_ss.faa"
        SeqIO.write(to_write, path, "fasta")
        print("Proteins written to: "+path)
    except:
        print("Error writing proteins to file!")
                          
def main():
    ### Logo ###
    if logo == True:
        print_logo()
        
    ### Run Scan ###
    # Some failsafes
    if seqs is None:
        print('Please provide at least one sequence with -s')
        sys.exit(1)
    cpus = str(threads)
    if threads == 0:
        cpus = str(os.cpu_count())
    if threads > os.cpu_count():
        print('Number of threads exceeds available CPUs, ' \
          'defaulting to maximum available CPUs')
        cpus = str(os.cpu_count())
    print('Using '+cpus+' threads...')
    if out != 'stdout':
        df = pd.DataFrame()
    
    # Loop over each input file
    for seq in seqs:      
        # Check if input is empty
        if os.path.isfile(seq) == False:
            print("No such file: "+seq)
            print()
            continue

        if os.path.getsize(seq) <= 10:
            print(seq+" is empty, skipping...")
            print()
            continue
        
        # Check molecule type and restrictions
        try:
            test = next(SeqIO.parse(seq,"fasta"))
            pass
        except:
            print(seq+" is not a fasta file, skipping...")
            print()
            continue
        
        if "E" in test._seq._data:
            molecule = 'aa'
            print(os.path.splitext(os.path.basename(seq))[0]+
                  ': protein input detected')
        else:
            molecule = 'dna'
            print(os.path.splitext(os.path.basename(seq))[0]+
                  ": DNA input detected")
        
        # If DNA input, get proteins
        quality = 'single'
        if molecule == 'dna':
            small_contig_list = []
            for r in SeqIO.parse(seq, format = "fasta"):
                if len(r.seq) < 10000:
                    small_contig_list.append(r.id)
            if len(small_contig_list) > 50:
                print(str(len(small_contig_list))+' contigs are < 10kbp')
                quality = 'meta'
            
            if len([1 for line in open(seq) if line.startswith(">")]) >= 1000:
                print(os.path.splitext(os.path.basename(seq))[0]+' has >= 1000 contigs')
                quality = 'meta'
            
            if lowq == True:
                quality = 'meta'
            
            if quality == 'meta':
                print('Switching Prodigal to anonymous mode')
                
            proteins = run_prodigal(seq, quality)
            if len(proteins) == 0:
                print("No CDS found")
                print()
                continue
        else: proteins = seq
        
        # Filter by domain, also shrinks fasta in memory                 
        proteins = run_hmmsearch(molecule, proteins, cpus, quality)
        if len(proteins) == 0:
            print("No TonB-dependent receptors found")
            print()
            continue
        hits = run_hmmscan(proteins, cpus, molecule)
        if hits.empty == True:
            continue
        
        if genloc == True:
            if molecule == 'aa':
                print("Cannot screen plasmids in protein fasta")
            else:
                hits = plasmid_screen(seq, cpus, hits)
                hits = mge_screen(seq, cpus, hits)
                    
        if fur == True:
            if molecule == 'aa':
                print("Cannot screen TFBS in protein fasta")
            else:
                hits = run_mast(hits, seq)
        
        if len(write) > 0:
            grab_proteins(proteins, hits, write, seq)
        
        # Tidying up
        hits = hits.drop(columns=['contig'])
        hits.insert(0, 'sample', os.path.splitext(os.path.basename(seq))[0])
        hits = hits.fillna(value = '-')
        hits = hits.set_index('sample') 
         
        if out == 'stdout':
            display(hits)
            print()
        else:
            df = df.append(hits)

    if out != 'stdout':
        df = df.fillna(value = '-')
        df.to_csv(out)
        print("Done, written to: "+out)

if __name__ == "__main__":
    main()