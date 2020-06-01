#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'sideroscanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import argparse
import os
import re
import sys
import textwrap
import warnings
from shutil import which, get_terminal_size
from subprocess import Popen, PIPE, DEVNULL
from argparse import RawTextHelpFormatter
from datetime import datetime
from io import StringIO
import pandas as pd
from Bio import SeqIO, SearchIO, BiopythonWarning
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from IPython.display import display

warnings.simplefilter('ignore', BiopythonWarning)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
        SideroScannner: A tool for annotating IROMPs in bacteria.

        By Tom Stanton - Schneiders Lab - University of Edinburgh
        For queries/feedback/comments/collaborations/chat/coffee donation:
        Email: T.D.Stanton@sms.ed.ac.uk --- Github: https://github.com/tomdstanton
        
        Please note:
            Curating the IROMP HMM library is an ongoing process.'''))

    parser.add_argument('-t', '--threads', type=int, default=os.cpu_count(),
                        help='Number of threads to use, (default: '+str(os.cpu_count())+')')

    parser.add_argument('-i', '--input', nargs='*', type=str,  required=True,
                        help='Path to fasta input, accepts multiple DNA or AA files.')

    parser.add_argument('-o', '--out', nargs='?', type=str,
                        const='SideroScanner_' + datetime.now().strftime("%d%m%y_%H%M") + '.csv',
                        help='''Save the output as a comma-separated file instead of STDOUT
Output path/file - (default: SideroScanner_DDMMYY_hhmm.csv)''')

    parser.add_argument('-g', '--genloc', nargs='?', type=str, const='90',
                        help='''If using contigs/assemblies, screen the input file
for plasmids/MGEs to determine genomic location of hits
If needed, enter value for blastn percid (default: 90)''')

    parser.add_argument('-f', '--flanks', nargs='?', type=int, const=3,
                        help='''Screen hit flanking CDS for virulence factors
Specifiy number of flanking CDS to screen - (default: 3)''')

    parser.add_argument('--tfbs', nargs='?', type=str, const='pwm',
                        help='''If using contigs/assemblies, scan the promoter regions of hits for TFBS using MAST
Provide a custom PSSM in meme format - (default: Fur pwm)''')

    parser.add_argument('-e', '--export', nargs='?', type=str,
                        const='seq',
                        help='''Write proteins from analysis to file
Exported protein path/file - (default: infile_SideroScanner.faa)''')

    parser.add_argument('--lowqual', nargs='?', type=str, const='',
                        help='''Adjusts gene prediction mode and filters only with the plug domain 
Provide a list of filenames this applies to  - (default: all inputs)''')

    parser.add_argument('--lib', type=str, default='iromps.hmm',
                        help='''Path to IROMP HMM library
Must be "pressed" in hmmer3 format - (default: iromps.hmm)''')

    parser.add_argument('--dbpath', type=str, default=sys.path[0],
                        help='''Path to blast databases for plasmid, MGE and flanking VF screens
(default: '''+os.getcwd()+')')

    if len(sys.argv)==1:
        parser.print_help()
        print('\n'+'Please provide at least one argument or -h for help'+'\n')
        sys.exit(1)
    return parser.parse_args()

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None

# def prepare_dbs():
#     def download(url):
#         out_file_path = url.split("/")[-1][:-3]
#         print('Downloading SEED Database from: {}'.format(url))
#         response = urllib2.urlopen(url)
#         compressed_file = StringIO.StringIO(response.read())
#         decompressed_file = gzip.GzipFile(fileobj=compressed_file)
#         with open(out_file_path, 'w') as outfile:
#             outfile.write(decompressed_file.read())
#         return
#     vfdb = download('http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz')
#
#
#     ices = http://202.120.12.136/ICEberg2/download/ICE_seq_all.fas
#     t4ss = http://202.120.12.136/ICEberg2/download/T4SS-type_ICE_seq_all.fas
#     aices = http://202.120.12.136/ICEberg2/download/AICE_seq_all.fas
#     imes = http://202.120.12.136/ICEberg2/download/IME_seq_all.fas
#     cimes = http://202.120.12.136/ICEberg2/download/CIME_seq_all.fas
#     plsdb = https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip
#     ["makeblastdb -in VFDB_setB_pro.fas -dbtype 'prot' -out vfdb"]
#
#     handle1 = fetch("http://www.uniprot.org/uniprot/", params).decode('utf-8').splitlines()
#     records = list(SeqIO.parse(handle1, "fasta"))
#     for r in records:
#         r.id = iromp
#         r.description = r.description.replace(r.name,'').lstrip()
#         sequences.append(r)
#     import gzip
    


def run_prodigal(infile, quality):
    print("Extracting proteins...")
    cmd = ['prodigal', '-i', infile, '-o', '/dev/null', '-a',
           '/dev/stdout', '-q', '-p', quality]
    prodigal = Popen(cmd, stdout=PIPE)
    return prodigal.communicate()[0].decode('utf-8')


def run_hmmsearch(molecule, infile, cpus, qual):
    if molecule == 'aa':
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus,
               parse_args().dbpath+'/PF07715.hmm', infile]
        hmmsearch = Popen(cmd, stdout=PIPE)
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        qresult = next(SearchIO.parse(hmmer_out, 'hmmer3-text'))
        records = list(SeqIO.parse(infile, "fasta"))
    else:
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus,
               parse_args().dbpath+'/PF07715.hmm', '-']
        hmmsearch = Popen(cmd, stdout=PIPE, stdin=PIPE)
        hmmsearch.stdin.write(infile.encode())
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        qresult = next(SearchIO.parse(hmmer_out, 'hmmer3-text'))
        records = SeqIO.parse(StringIO(infile), 'fasta')
    filter1 = ''
    for r in records:
        if r.id in qresult.hit_keys:
            filter1 = filter1 + r.format("fasta")
    print('Filtered '+str(filter1.count('>'))+' proteins with PF07715.hmm')
    if qual == 'meta' or len(filter1) == 0:
        return filter1
    else:
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus,
               parse_args().dbpath+'/PF00593.hmm', '-']
        hmmsearch = Popen(cmd, stdin=PIPE,
                          stdout=PIPE)
        hmmsearch.stdin.write(filter1.encode())
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        qresult = next(SearchIO.parse(hmmer_out, 'hmmer3-text'))
        records = SeqIO.parse(StringIO(filter1), 'fasta')
        filter2 = ''
        for r in records:
            if r.id in qresult.hit_keys:
                filter2 = filter2 + r.format("fasta")
        print('Filtered '+str(filter2.count('>'))+' proteins with PF00593.hmm')
        return filter2


def run_hmmscan(infile, cpus, molecule):
    print("Scanning against HMM library...")
    hmmscan_df = pd.DataFrame(columns=['contig', 'query', 'hit', 'hit_range',
                                       'score', 'evalue'])
    cmd = ['hmmscan', '-E', '1e-50', '--cpu', cpus, parse_args().lib, '-']
    hmmscan = Popen(cmd, stdin=PIPE, stdout=PIPE)
    hmmscan.stdin.write(infile.encode())
    hmmer_out = StringIO(hmmscan.communicate()[0].decode('utf-8'))
    for qresult in SearchIO.parse(hmmer_out, "hmmer3-text"):
        if len(qresult.hits) > 0:
            hmmscan_df = hmmscan_df.append(
                {'contig': qresult._id.rsplit('_', 1)[0], 'query': qresult._id,
                 'hit': qresult.hits[0]._id, 'score': qresult.hits[0].bitscore,
                 'evalue': "%.3g" % qresult.hits[0].evalue,
                 'hit_range': str(qresult.hits[0].hsps[0].hit_range[0]) \
                              + '-' + str(qresult.hits[0].hsps[0].hit_range[1])},
                ignore_index=True)
    if hmmscan_df.empty == True:
        print("No proteins annotated")
        return hmmscan_df
    weight_df = pd.DataFrame(columns=['query', 'mass(kDa)'])
    for r in SeqIO.parse(StringIO(infile), 'fasta'):
        X = ProteinAnalysis((r._seq._data.replace('*', '')))
        weight_df = weight_df.append({'query': r.id, 'mass(kDa)':round((X.molecular_weight()/1000),1)}, ignore_index=True)
    hmmscan_df = hmmscan_df.merge(weight_df, on='query', how='left')
    if molecule == 'aa':
        return hmmscan_df
    else:
        bed_df = pd.DataFrame(columns=['query', 'start', 'end', 'strand'])
        for r in SeqIO.parse(StringIO(infile), 'fasta'):
            bed = re.split('\#|\s', (r.description + "\n").replace(' ', ''))[0:4]
            bed_df = bed_df.append(pd.Series(
                bed, index=['query', 'start', 'end', 'strand']),
                ignore_index=True)
        return pd.merge(hmmscan_df, bed_df, on='query')


def plasmid_screen(infile, cpus, hits):
    print("Screening for plasmids...")
    cmd = ['blastn',
           '-query', infile,
           '-db', parse_args().dbpath+'/plsdb.fna',
           '-culling_limit', '1',
           '-outfmt', '6',
           '-max_hsps', '1',
           '-evalue', '1e-50',
           '-num_threads', cpus,
           '-perc_identity', parse_args().genloc,
           '-qcov_hsp_perc', '10',
           '-max_target_seqs', '1']
    blast = Popen(cmd, stdout=PIPE, stderr=DEVNULL)
    blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    plasmid_df = pd.DataFrame(columns=['contig', 'plasmid/mge', 'span'])
    for qresult in SearchIO.parse(blast_out, "blast-tab"):
        for h in qresult.hits:
            plasmid_df = plasmid_df.append(
                {'contig': h.query_id, 'plasmid/mge': h.id,
                 'span': str(h.hsps[0].query_range[0]) + '-' + str(h.hsps[0].query_range[1])}, ignore_index=True)
    hits = hits.merge(plasmid_df, how='left')
    if plasmid_df.empty == True:
        print("No plasmids found")
    return hits.fillna(value={'plasmid/mge': 'Chromosome'})


def mge_screen(infile, cpus, plasmids):
    print("Screening for MGEs...")
    cmd = ['blastn',
           '-query', infile,
           '-db', parse_args().dbpath+'/mge',
           '-culling_limit', '1',
           '-outfmt', '6',
           '-max_hsps', '1',
           '-evalue', '1e-50',
           '-num_threads', cpus,
           '-perc_identity', parse_args().genloc]
    blast = Popen(cmd, stdout=PIPE)
    blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    mge_df = pd.DataFrame(columns=['contig', 'plasmid/mge', 'span'])
    for qresult in SearchIO.parse(blast_out, "blast-tab"):
        for h in qresult.hits:
            name = h.id.split('|')[2]
            mge_df = mge_df.append(
                {'contig': h.query_id, 'plasmid/mge': name,
                 'span': h.hsps[0].query_range}, ignore_index=True)
    chromosome_contigs = list(plasmids.loc[plasmids['plasmid/mge'] == 'Chromosome', 'contig'])
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
                                range(mge_dict[key][0], mge_dict[key][1])) == True:
                    plasmids.at[row.Index, 'plasmid/mge'] = key
                    plasmids.at[row.Index, 'span'] = \
                        str(mge_dict[key][0]) + '-' + str(mge_dict[key][1])
    return plasmids


def flank_screen(infile, cpus, hits, cds):
    flanking_cds = ''
    for r in SeqIO.parse(StringIO(infile), 'fasta'):
        if r.id in hits['query'].values:
            for i in range(cds):
                upstream_cds = r.id.rsplit('_', 1)[0]+'_'+str(int(r.id.rsplit('_', 1)[1])-int(i+1))
                downstream_cds = r.id.rsplit('_', 1)[0]+'_'+str(int(r.id.rsplit('_', 1)[1])+int(i+1))
                for rec in SeqIO.parse(StringIO(infile), 'fasta'):
                    if rec.id == upstream_cds:
                        flanking_cds = flanking_cds + '>' + str(r.id) + ' upstream '+ str(i+1) + '\n' + str(rec.seq) + '\n'
                    if rec.id == downstream_cds:
                        flanking_cds = flanking_cds + '>' + str(r.id) + ' downstream '+ str(i+1) + '\n' + str(rec.seq) + '\n'
    print('Screening '+str(cds)+' upstream and downstream CDS...')
    cmd = ['blastp',
           '-db', parse_args().dbpath+'/flankdb',
           '-outfmt', '5',
           '-max_hsps', '1',
           '-evalue', '1e-130',
           '-num_threads', cpus,
           '-query', '-']
    blast = Popen(cmd, stdout=PIPE, stdin=PIPE)
    blast.stdin.write(flanking_cds.encode())
    blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    u_df = pd.DataFrame(columns=['query', 'upstream'])
    d_df = pd.DataFrame(columns=['query', 'downstream'])
    for qresult in SearchIO.parse(blast_out, "blast-xml"):
        if len(qresult.hits) > 0:
            gene_name = qresult.hsps[0].hit_description.split('(', 1)[1].split(')')[0]
            if 'upstream' in qresult.hsps[0].query_description:
                u_df = u_df.append({'query': qresult.hsps[0].query_id,'upstream':gene_name},ignore_index=True)
            if 'downstream' in qresult.hsps[0].query_description:
                d_df = d_df.append({'query': qresult.hsps[0].query_id,'downstream':gene_name},ignore_index=True) 
    u_df = u_df.groupby(['query'])['upstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    d_df = d_df.groupby(['query'])['downstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    hits = hits.merge(u_df, on='query', how='left')
    hits = hits.merge(d_df, on='query', how='left')
    return hits


def run_mast(hits, pwm, infile):
    print('Searching for motif hits in promoter regions...')
    hits.loc[hits['strand'] == '1', 'flank_start'] = hits['start'].astype(int) - 450
    hits.loc[hits['strand'] == '1', 'flank_end'] = hits['start'].astype(int) + 150
    hits.loc[hits['strand'] == '-1', 'flank_start'] = hits['end'].astype(int) - 150
    hits.loc[hits['strand'] == '-1', 'flank_end'] = hits['end'].astype(int) + 450
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
        promoter_reg = promoter_reg + '>' + query + ':' + str(f_start) + '-' + str(f_end) + '\n' + str(flank) + '\n'
    cmd = ['mast', '-hit_list', pwm, '-']
    mast = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=DEVNULL)
    mast.stdin.write(promoter_reg.encode())
    mast_out = StringIO(mast.communicate()[0].decode('utf-8'))
    mast_df = pd.read_csv(mast_out, sep=" ", header=1)
    if mast_df.empty == True:
        print("No TFBS found")
        return hits.drop(columns=['flank_start', 'flank_end'])
    else:
        mast_df = mast_df[:-1]
        mast_df[['query', 'flank']] = mast_df['#'].str.split(":", expand=True)
        mast_df[['flank_start', 'flank_end']] = mast_df['flank'].str.split("-", expand=True)
        mast_df = mast_df.drop(columns=['#', 'id', 'score', 'flank', 'hit_end', 'sequence_name'])
        mast_df = mast_df.rename(columns={"(strand+/-)motif": "motif",
                                          "alt_id": "motif_start",
                                          "hit_start": "motif_end",
                                          "hit_p-value": "motif_p-value"})
        mast_df['motif_start'] = (mast_df['motif_start'].astype(int) + mast_df['flank_start'].astype(int)).astype(str)
        mast_df['motif_end'] = (mast_df['motif_end'].astype(int) + mast_df['flank_start'].astype(int)).astype(str)
        mast_df = mast_df.drop(columns=['flank_start', 'flank_end'])
        hits = hits.drop(columns=['flank_start', 'flank_end'])
        hits = hits.merge(mast_df, how='left', on='query')
        return hits
    

def grab_proteins(infile, hits, seq):
    to_write = []
    for r in SeqIO.parse(StringIO(infile), 'fasta'):
        if r.id in hits['query'].values:
            write_df = hits[hits['query'].str.match(r.id)]
            r.id = write_df['query'].values[0] + '_' + write_df['hit'].values[0]
            r.description = ''
            to_write.append(r)
    if parse_args().export == 'seq':
        filename = os.path.splitext(seq)[0]+'_SideroScanner.faa'
    else:
        filename = parse_args().export
    SeqIO.write(to_write, filename, "fasta")
    print("Proteins written to: " + filename)


def main():
    print('-'*int(get_terminal_size()[0]))
    for d in ['hmmscan', 'hmmsearch', 'blastp', 'blastn', 'prodigal', 'mast']:
        if is_tool(d) == False:
            print('ERROR: '+d+' not found')
            sys.exit(1)
    
    if parse_args().input is None:
        print('Please provide at least one sequence with -s')
        sys.exit(1)

    print(__title__+': '+__version__)
    print('Your system is '+os.uname()[0])
    print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))

    if parse_args().threads > os.cpu_count():
        threads = str(os.cpu_count())
        print('Number of threads exceeds available CPUs, will use: ' + threads)
    else:
        threads = str(parse_args().threads)
        print('Using ' + threads + ' threads...')

    if parse_args().out is not None:
        output_df = pd.DataFrame()

    # Loop over each input file
    for seq in parse_args().input:
        print('-'*int(get_terminal_size()[0]))
        if not os.path.isfile(seq):
            print("No such file: " + seq)
            continue

        if os.path.getsize(seq) <= 10:
            print(seq + " is too small, skipping...")
            continue

        # Check molecule type and restrictions
        try:
            test = next(SeqIO.parse(seq, "fasta"))
            pass
        except:
            print(seq + " is not a fasta file, skipping...")
            continue
        if "E" in test._seq._data:
            molecule = 'aa'
            print('Protein input detected: '+os.path.splitext(os.path.basename(seq))[0])
        else:
            molecule = 'dna'
            print('DNA input detected: '+os.path.splitext(os.path.basename(seq))[0])

        # If DNA input, get proteins
        quality = 'single'
        if molecule == 'dna':
            small_contig_list = []
            for r in SeqIO.parse(seq, format="fasta"):
                if len(r.seq) < 10000:
                    small_contig_list.append(r.id)

            if len(small_contig_list) > 50:
                print(str(len(small_contig_list)) + ' contigs are < 10kbp')
                quality = 'meta'

            if len([1 for line in open(seq) if line.startswith(">")]) >= 1000:
                print(os.path.splitext(os.path.basename(seq))[0] + ' has >= 1000 contigs')
                quality = 'meta'

            if parse_args().lowqual is not None:
                if len(parse_args().lowqual) == 0:
                    quality = 'meta'
                if seq in parse_args().lowqual:
                    quality = 'meta'

            if quality == 'meta':
                print('Switching Prodigal to anonymous mode')

            proteins = run_prodigal(seq, quality)
            if len(proteins) == 0:
                print("No CDS found")
                continue
        else:
            proteins = seq
             
        iromps = run_hmmsearch(molecule, proteins, threads, quality)
        if len(iromps) == 0:
            print('No TonB-dependent receptors found')
            continue
        
        hits = run_hmmscan(iromps, threads, molecule)
        if hits.empty is True:
            continue
        
        if parse_args().flanks is not None:
            if molecule == 'aa':
                print('Cannot screen flanking CDS in protein fasta')
            else:
                hits = flank_screen(proteins, threads, hits, parse_args().flanks)

        del proteins

        if parse_args().genloc is not None:
            if molecule == 'aa':
                print('Cannot screen plasmids in protein fasta')
            else:
                hits = plasmid_screen(seq, threads, hits)
                hits = mge_screen(seq, threads, hits)

        if parse_args().tfbs is not None:
            if molecule == 'aa':
                print("Cannot screen promoter regions in protein fasta")
            else:
                hits = run_mast(hits, parse_args().tfbs, seq)

        if parse_args().export is not None:
            grab_proteins(iromps, hits, seq)

        # Tidying up
        hits = hits.drop(columns=['contig'])
        hits.insert(0, 'sample', os.path.splitext(os.path.basename(seq))[0])
        hits = hits.fillna(value='-')
        hits = hits.set_index('sample')

        if parse_args().out is None:
            display(hits)
        else:
            output_df = output_df.append(hits)

    if parse_args().out is not None:
        output_df = output_df.fillna(value='-')
        output_df.to_csv(parse_args().out)
        print("Done, written to: " + parse_args().out)


if __name__ == "__main__":
    main()
