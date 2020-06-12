#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'sideroscanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import argparse
import os
import re
import sys
import warnings
from shutil import which, get_terminal_size
from subprocess import Popen, PIPE, DEVNULL, run
from argparse import RawTextHelpFormatter
from datetime import datetime
from io import StringIO
import pandas as pd
from Bio import SeqIO, SearchIO, BiopythonWarning
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from IPython.display import display
from build_dbs import fetch

warnings.simplefilter('ignore', BiopythonWarning)
pathname = os.path.dirname(sys.argv[0])
full_path = os.path.abspath(pathname)

def parse_args():
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=RawTextHelpFormatter,
                                     usage="./sideroscanner.py -i [-f -l -o -e]",
                                     description='''
        sideroscannner: a tool for annotating IROMPs in bacteria
        ========================================================''')
    group = parser.add_argument_group("Options")

    group.add_argument('-i', metavar='-', nargs='*', type=str,
                        help='''| path/to/(i)nput/fasta
-----------------------------------------------''')
    group.add_argument('-o', metavar='-', nargs='?', type=str,
                        const='sideroscanner_' + datetime.now().strftime("%d%m%y_%H%M") + '.csv',
                        help='''| (o)utput file.csv instead of STDOUT
    [optional: path/to/(o)utput/file]
    [default: sideroscanner_DDMMYY_hhmm.csv]
-----------------------------------------------''')
    group.add_argument('-l', metavar='int', nargs='?', type=str, const='90',
                        help='''| determine genomic (l)ocation of hits
    [optional: blastn percid]
    [default: 90]
-----------------------------------------------''')
    group.add_argument('-f', metavar='int', nargs='?', type=int, const=3,
                        help='''| (f)lanking CDS screen
    [optional: number of up/downstream CDS]
    [default: 3]
-----------------------------------------------''')
    group.add_argument('-e', metavar='-', nargs='?', type=str,
                        const='seq',
                        help='''| (e)xport annotated proteins
    [optional: path/to/export/fasta]
    [default: *infile*_sideroscanner.faa]
-----------------------------------------------''')
    group.add_argument('-t', metavar='int', type=int, default=os.cpu_count(),
                       help='''| number of (t)hreads to use
    [default: max available]
-----------------------------------------------''')
    group.add_argument('--lowqual', metavar='-', nargs='?', type=str, const='',
                        help='''| (1) 'meta' CDS prediction AND (2) filters with plug domain only
    [optional: path/to/(low)/(qual)ity/input/fasta]
    [default: all inputs]
-----------------------------------------------''')
    group.add_argument('--lib', metavar='hmm', type=str, default=full_path + '/databases/iromps.hmm',
                        help='''| path/to/custom/HMM
    [default: ''' + full_path + '''/databases/iromps.hmm]
-----------------------------------------------''')
    group.add_argument('--dbpath', metavar='path', type=str, default=full_path + '/databases/',
                        help='''| path/to/db/
    [default: ''' + full_path + '''/databases/]
-----------------------------------------------''')
    group.add_argument('-v', action='store_true',
                       help='''| show version and exit
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''| show this help message and exit''')

    if len(sys.argv)==1:
        parser.print_help()
        print('''
Please provide at least one argument or -h for help
        ''')
        sys.exit(1)
    return parser.parse_args()


dbpath = parse_args().dbpath
lib = parse_args().lib

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None

def run_prodigal(infile, quality):
    print("Extracting proteins...")
    cmd = ['prodigal', '-i', infile, '-o', '/dev/null', '-a',
           '/dev/stdout', '-q', '-p', quality]
    prodigal = Popen(cmd, stdout=PIPE)
    return prodigal.communicate()[0].decode('utf-8')


def run_hmmsearch(molecule, infile, cpus, qual):
    if molecule == 'aa':
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus,
               dbpath + '/PF07715.hmm', infile]
        hmmsearch = Popen(cmd, stdout=PIPE)
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        q = next(SearchIO.parse(hmmer_out, 'hmmer3-text'))
        records = list(SeqIO.parse(infile, "fasta"))
    else:
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus,
               dbpath + '/PF07715.hmm', '-']
        hmmsearch = Popen(cmd, stdout=PIPE, stdin=PIPE)
        hmmsearch.stdin.write(infile.encode())
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        q = next(SearchIO.parse(hmmer_out, 'hmmer3-text'))
        records = SeqIO.parse(StringIO(infile), 'fasta')
    filter1 = ''
    for r in records:
        if r.id in q.hit_keys:
            filter1 = filter1 + r.format("fasta")
    print('Filtered ' + str(filter1.count('>')) + ' proteins with PF07715.hmm')
    if qual == 'meta' or len(filter1) == 0:
        return filter1
    else:
        cmd = ['hmmsearch', '--cut_tc', '--cpu', cpus,
               dbpath + '/PF00593.hmm', '-']
        hmmsearch = Popen(cmd, stdin=PIPE,
                          stdout=PIPE)
        hmmsearch.stdin.write(filter1.encode())
        hmmer_out = StringIO(hmmsearch.communicate()[0].decode('utf-8'))
        q = next(SearchIO.parse(hmmer_out, 'hmmer3-text'))
        records = SeqIO.parse(StringIO(filter1), 'fasta')
        filter2 = ''
        for r in records:
            if r.id in q.hit_keys:
                filter2 = filter2 + r.format("fasta")
        print('Filtered ' + str(filter2.count('>')) + ' proteins with PF00593.hmm')
        return filter2


def run_hmmscan(infile, cpus, molecule):
    print("Scanning against HMM library...")
    hmmscan_df = pd.DataFrame(columns=['contig', 'query', 'hit', #'hit_range',
                                       'score', 'evalue'])
    cmd = ['hmmscan', '-E', '1e-50', '--cpu', cpus, lib, '-']
    hmmscan = Popen(cmd, stdin=PIPE, stdout=PIPE)
    hmmscan.stdin.write(infile.encode())
    hmmer_out = StringIO(hmmscan.communicate()[0].decode('utf-8'))
    for q in SearchIO.parse(hmmer_out, "hmmer3-text"):
        if len(q.hits) > 0:
            hmmscan_df = hmmscan_df.append(
                {'contig': q._id.rsplit('_', 1)[0], 'query': q._id,
                 'hit': q.hits[0]._id, 'score': q.hits[0].bitscore,
                 'evalue': "%.3g" % q.hits[0].evalue #,
                 #'hit_range': str(q.hits[0].hsps[0].hit_range[0]) \
                 #             + '-' + str(q.hits[0].hsps[0].hit_range[1])
                 },
                ignore_index=True)
    if hmmscan_df.empty == True:
        print("No proteins annotated")
        return hmmscan_df
    len_mass_df = pd.DataFrame(columns=['query','len', 'mass(kDa)'])
    for r in SeqIO.parse(StringIO(infile), 'fasta'):
        L = str(len(r.seq))
        if molecule == 'dna':
            if '00' not in re.search(r'partial=(.*);start_type', r.description).group():
                L = L+'-partial'
        X = ProteinAnalysis((r._seq._data.replace('*', '')))
        len_mass_df = len_mass_df.append({'query': r.id,
                                          'len': L,
                                          'mass(kDa)': round((X.molecular_weight() / 1000), 1)},
                                     ignore_index=True)
    hmmscan_df = hmmscan_df.merge(len_mass_df, on='query', how='left')
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
           '-db', dbpath+'plsdb.fna',
           '-culling_limit', '1',
           '-outfmt', '6',
           '-max_hsps', '1',
           '-evalue', '1e-50',
           '-num_threads', cpus,
           '-perc_identity', '90',
           '-qcov_hsp_perc', '10',
           '-max_target_seqs', '1']
    blast = Popen(cmd, stdout=PIPE, stderr=DEVNULL)
    blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    plasmid_df = pd.DataFrame(columns=['contig', 'plasmid/mge', 'span'])
    for q in SearchIO.parse(blast_out, "blast-tab"):
        for h in q.hits:
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
           '-db', dbpath+'mge',
           '-culling_limit', '1',
           '-outfmt', '6',
           '-max_hsps', '1',
           '-evalue', '1e-50',
           '-num_threads', cpus,
           '-perc_identity', parse_args().l]
    blast = Popen(cmd, stdout=PIPE)
    blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
    mge_df = pd.DataFrame(columns=['contig', 'plasmid/mge', 'span'])
    for q in SearchIO.parse(blast_out, "blast-tab"):
        for h in q.hits:
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
                upstream_cds = r.id.rsplit('_', 1)[0] + '_' + str(int(r.id.rsplit('_', 1)[1]) - int(i + 1))
                downstream_cds = r.id.rsplit('_', 1)[0] + '_' + str(int(r.id.rsplit('_', 1)[1]) + int(i + 1))
                for rec in SeqIO.parse(StringIO(infile), 'fasta'):
                    if rec.id == upstream_cds:
                        flanking_cds = flanking_cds + '>' + str(r.id) + ' upstream ' + str(i + 1) + '\n' + str(
                            rec.seq) + '\n'
                    if rec.id == downstream_cds:
                        flanking_cds = flanking_cds + '>' + str(r.id) + ' downstream ' + str(i + 1) + '\n' + str(
                            rec.seq) + '\n'
    print('Screening ' + str(cds) + ' upstream and downstream CDS...')
    cmd = ['blastp',
           '-db', dbpath+'flankdb',
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
    for q in SearchIO.parse(blast_out, "blast-xml"):
        if len(q.hits) > 0:
            gene_name = re.sub("\s*\[[^[]*\]$", '', q.hsps[0].hit_description).strip()
            if gene_name.startswith('('):
                gene_name = gene_name.split('(', 1)[1].split(')')[0]
            if 'upstream' in q.hsps[0].query_description:
                u_df = u_df.append({'query': q.hsps[0].query_id, 'upstream': gene_name}, ignore_index=True)
            if 'downstream' in q.hsps[0].query_description:
                d_df = d_df.append({'query': q.hsps[0].query_id, 'downstream': gene_name}, ignore_index=True)
    u_df = u_df.groupby(['query'])['upstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    d_df = d_df.groupby(['query'])['downstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    hits = hits.merge(u_df, on='query', how='left')
    hits = hits.merge(d_df, on='query', how='left')
    return hits


# def run_mast(hits, pwm, infile):
#     print('Searching for motif hits in promoter regions...')
#     hits.loc[hits['strand'] == '1', 'flank_start'] = hits['start'].astype(int) - 450
#     hits.loc[hits['strand'] == '1', 'flank_end'] = hits['start'].astype(int) + 150
#     hits.loc[hits['strand'] == '-1', 'flank_start'] = hits['end'].astype(int) - 150
#     hits.loc[hits['strand'] == '-1', 'flank_end'] = hits['end'].astype(int) + 450
#     chrs = {}
#     for r in SeqIO.parse(infile, "fasta"):
#         chrs[r.id] = r.seq
#     promoter_reg = ''
#     for index, row in hits.iterrows():
#         query = row['query']
#         chracc = row['contig']
#         f_start = int(row['flank_start'])
#         f_end = int(row['flank_end'])
#         flank = chrs[chracc][f_start:f_end]
#         promoter_reg = promoter_reg + '>' + query + ':' + str(f_start) + '-' + str(f_end) + '\n' + str(flank) + '\n'
#     cmd = ['mast', '-hit_list', pwm, '-']
#     mast = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=DEVNULL)
#     mast.stdin.write(promoter_reg.encode())
#     mast_out = StringIO(mast.communicate()[0].decode('utf-8'))
#     mast_df = pd.read_csv(mast_out, sep=" ", header=1)
#     if mast_df.empty == True:
#         print("No TFBS found")
#         return hits.drop(columns=['flank_start', 'flank_end'])
#     else:
#         mast_df = mast_df[:-1]
#         mast_df[['query', 'flank']] = mast_df['#'].str.split(":", expand=True)
#         mast_df[['flank_start', 'flank_end']] = mast_df['flank'].str.split("-", expand=True)
#         mast_df = mast_df.drop(columns=['#', 'id', 'score', 'flank', 'hit_end', 'sequence_name'])
#         mast_df = mast_df.rename(columns={"(strand+/-)motif": "motif",
#                                           "alt_id": "motif_start",
#                                           "hit_start": "motif_end",
#                                           "hit_p-value": "motif_p-value"})
#         mast_df['motif_start'] = (mast_df['motif_start'].astype(int) + mast_df['flank_start'].astype(int)).astype(str)
#         mast_df['motif_end'] = (mast_df['motif_end'].astype(int) + mast_df['flank_start'].astype(int)).astype(str)
#         mast_df = mast_df.drop(columns=['flank_start', 'flank_end'])
#         hits = hits.drop(columns=['flank_start', 'flank_end'])
#         hits = hits.merge(mast_df, how='left', on='query')
#         return hits


def grab_proteins(infile, hits, seq):
    to_write = []
    for r in SeqIO.parse(StringIO(infile), 'fasta'):
        if r.id in hits['query'].values:
            write_df = hits[hits['query'].str.match(r.id)]
            r.id = write_df['query'].values[0] + '_' + write_df['hit'].values[0]
            r.description = ''
            to_write.append(r)
    if parse_args().e == 'seq':
        filename = os.path.splitext(seq)[0] + '_sideroscanner.faa'
    else:
        filename = parse_args().e
    SeqIO.write(to_write, filename, "fasta")
    print("Proteins written to: " + filename)


def main():
    # Some quick sanity checks
    # Check if tools are installed
    for d in ['hmmscan', 'hmmsearch', 'blastp', 'blastn', 'prodigal']:
        if is_tool(d) == False:
            print('ERROR: ' + d + ' not found')
            sys.exit(1)

    # Check database folder
    if not os.path.isdir(dbpath):
        print(dbpath + 'is not an existing directory')
        sys.exit(1)
    else:
        plug = dbpath + 'PF07715.hmm'
        if not os.path.isfile(plug):
            print (plug+' not found, fetching...')
            with open(plug, mode='wb') as localfile:
                localfile.write(fetch('https://pfam.xfam.org/family/PF07715/hmm'))
        receptor = dbpath + 'PF00593.hmm'
        if not os.path.isfile(receptor):
            print(receptor + ' not found, fetching...')
            with open(receptor, mode='wb') as localfile:
                localfile.write(fetch('https://pfam.xfam.org/family/PF00593/hmm'))

    # Check library HMM
    if not os.path.isfile(lib):
        print(lib + ' is not a valid file')
        sys.exit(1)
    else:
        for i in [".h3i", ".h3p", ".h3f", ".h3m"]:
            if not os.path.isfile(lib + i):
                print(lib + ' is not pressed, allow me...')
                run(["hmmpress", lib])
                break

    # Sanity checks passed, assuming good to go
    print('-' * int(get_terminal_size()[0]))
    print(__title__ + ': ' + __version__)
    print('Your system is ' + os.uname()[0])
    if 'Linux' not in os.uname()[0]:
        print('Warning: sideroscanner has not been tested on '+os.uname()[0])
    print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))

    if parse_args().t > os.cpu_count():
        threads = str(os.cpu_count())
        print('Number of threads exceeds available CPUs, will use: ' + threads)
    else:
        threads = str(parse_args().t)
        print('Using ' + threads + ' threads...')

    if parse_args().o is not None:
        out_df = pd.DataFrame()

    # Loop over each input file
    for seq in parse_args().i:
        print('-' * int(get_terminal_size()[0]))
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
            print('Protein input detected: ' + os.path.splitext(os.path.basename(seq))[0])
        else:
            molecule = 'dna'
            print('DNA input detected: ' + os.path.splitext(os.path.basename(seq))[0])

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

        if parse_args().f is not None:
            if molecule == 'aa':
                print('Cannot screen flanking CDS in protein fasta')
            else:
                hits = flank_screen(proteins, threads, hits, parse_args().f)

        del proteins

        if parse_args().l is not None:
            if molecule == 'aa':
                print('Cannot screen plasmids in protein fasta')
            else:
                hits = plasmid_screen(seq, threads, hits)
                hits = mge_screen(seq, threads, hits)

        # if parse_args().tfbs is not None:
        #     if molecule == 'aa':
        #         print("Cannot screen promoter regions in protein fasta")
        #     else:
        #         hits = run_mast(hits, parse_args().tfbs, seq)

        if parse_args().e is not None:
            grab_proteins(iromps, hits, seq)

        # Tidying up
        hits = hits.drop(columns=['contig'])
        hits.insert(0, 'sample', os.path.splitext(os.path.basename(seq))[0])
        hits = hits.fillna(value='-')
        hits = hits.set_index('sample')

        if parse_args().o is None:
            display(hits)
        else:
            out_df = out_df.append(hits)

    if parse_args().o is not None:
        out_df = out_df.fillna(value='-')
        out_df.to_csv(parse_args().o)
        print("Done, written to: " + parse_args().o)


if __name__ == "__main__":
    main()
