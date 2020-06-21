#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import sys, os, argparse, re
import pandas as pd
from argparse import RawTextHelpFormatter
from Bio import SeqIO, Entrez
from shutil import get_terminal_size
from datetime import datetime

from scripts.blast import run_blastp, run_makeblastdb
from scripts.cdhit import run_cdhit
from scripts.fetch import fetch
from scripts.hmmer3 import run_hmmbuild
from scripts.trimal import run_trimal
from scripts.muscle import run_muscle

pathname = os.path.dirname(sys.argv[0])
full_path = os.path.abspath(pathname)
cwd = os.getcwd()
threads = str(os.cpu_count())
sys.path.append(full_path)

def parse_args():
    parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter,
                                     usage="./sideroscanner.py -i genome.fna",
                                     description='''
            SideroScannner: a tool for annotating IROMPs in bacteria
            ========================================================''')
    group = parser.add_argument_group("Options")

    group.add_argument('-c', metavar='csv', type=str,
                       default=full_path + '/databases/irompdb/iromps.csv',
                       help='''| path/to/iromp/(c)sv
    [default: ''' + full_path + '''/databases/irompdb/iromps.csv]
-----------------------------------------------''')
    group.add_argument('--dbpath', metavar='path', type=str,
                       default=full_path + '/databases/',
                       help='''| path/to/db/
    [default: ''' + full_path + '''/databases/]
-----------------------------------------------''')
    group.add_argument('-w', metavar='int', type=int, default=2,
                       help='''| protein length (w)indow for blastp
    [default: 3]
-----------------------------------------------''')
    group.add_argument('-e', metavar='str', type=str, default='1e-250',
                       help='''| (e)value for blastp
    [default: 1e-250]
    -----------------------------------------------''')
    group.add_argument('-t', metavar='int', type=str, default=threads,
                       help='''| number of (t)hreads to use
    [default: ''' + threads + ''']
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''| show this help message and exit''')
    return parser.parse_args()


dbpath = parse_args().dbpath
window = parse_args().w
evalue = parse_args().e
iromp_csv = parse_args().c

if parse_args().t > threads:
    print('Number of threads exceeds available CPUs, will use: ' + threads)
else:
    threads = str(parse_args().t)


def get_fastadb():
    url = 'http://www.uniprot.org/uniprot/'
    params = {'query': 'taxonomy:"Bacteria [2]" AND '
                       'family:"tonb-dependent receptor family" AND '
                       'locations:(location:"Cell outer membrane [SL-0040]") '
                       'NOT partial NOT fragment', 'format': 'fasta'}
    print("Fetching TonB-depenent receptors")
    fetch(url, params, dbpath + 'irompdb/iromps.faa')
    return dbpath + 'irompdb/iromps.faa'


def seed_blast(in_file, blastdb, name, db_df):
    hit_acc = []
    for q in run_blastp(in_file, blastdb, evalue, threads):
        if len(q.hits) > 0:
            print('%i hits found' % len(q.hits))
            up = q.seq_len + window
            down = q.seq_len - window
            for h in q.hits:
                if down <= h.seq_len <= up:
                    hit_acc.append(h.blast_id)
            print("%i hits filtered between %i-%i aa" % (len(hit_acc), down, up))
            hit_df = db_df[db_df['accession'].isin(hit_acc)]
            hit_df.to_excel(writer, sheet_name=name, index=False)
            hits = ''
            for h in hit_acc:
                hits = hits + db_dict[h].format("fasta")
    return hits


def make_db_df(in_file):
    print('Creating DataFrame from ' + in_file + '...')
    data = []
    for r in SeqIO.parse(in_file, 'fasta'):
        if re.search(r'GN=(.*) PE=', r.description):
            gene = re.search(r'GN=(.*) PE=', r.description).group(1)
        else:
            gene = 'NA'
        data.append(r.id.split('|')[1] + ',' + gene + ',' +
                    re.search(r' (.*) OS=', r.description).group(1).replace(",", " ") + ',' +
                    re.search(r'OS=(.*) OX=', r.description).group(1).replace(",", " "))
    print('Done: %i total proteins in database' % len(data))
    return pd.DataFrame([sub.split(",") for sub in data],
                        columns=['accession', 'gene', 'description', 'organism'])


def fetch_seed(acc, name):
    print("Fetching " + name + " from NCBI...")
    Entrez.email = "tomdstanton@gmail.com"
    return Entrez.efetch(
        db='protein', id=acc, rettype="fasta",
        retmode="text").read().rstrip('\n')

print('-' * int(get_terminal_size()[0]))
print(__title__ + ': ' + __version__)
print('Your system is ' + os.uname()[0])
if 'Linux' not in os.uname()[0]:
    print('Warning: sideroscanner has not been tested on ' + os.uname()[0])
print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))
print('Using ' + threads + ' threads...')

if not os.path.isfile(dbpath + 'irompdb/iromps_nr.faa'):
    if not os.path.isfile(dbpath + 'irompdb/iromps.faa'):
        get_fastadb()
    else:
        fastadb = run_cdhit(dbpath + 'irompdb/iromps.faa', 1)
else:
    fastadb = dbpath + 'irompdb/iromps_nr.faa'

for i in [".psq", ".psi", ".psd", ".pog", ".pin", ".phr"]:
    if not os.path.isfile(dbpath+'irompdb/irompdb'+i):
        blastdb = run_makeblastdb(fastadb, 'prot', dbpath + 'irompdb/irompdb')
        break
    else:
        blastdb = dbpath + 'irompdb/irompdb'

db_df = make_db_df(fastadb)
db_dict = SeqIO.to_dict(SeqIO.parse(fastadb, 'fasta'))
seed_df = pd.read_csv(iromp_csv)
writer = pd.ExcelWriter(dbpath + 'irompdb/seed_alignment_blastp_' +
                        evalue + '_%i.xlsx' % window, engine='xlsxwriter')
for row in seed_df.itertuples():
    print('-' * int(get_terminal_size()[0]))
    iromp = fetch_seed(row.acc, row.protein)
    hits = seed_blast(iromp, blastdb, row.protein, db_df)
    print('-' * int(get_terminal_size()[0]))
    alignment = run_muscle(hits, dbpath + 'irompdb/' + row.protein + '_aligned.faa')
    trimmed = run_trimal(alignment)
    print('-' * int(get_terminal_size()[0]))
    hmm = run_hmmbuild(trimmed, row.protein, dbpath+'irompdb/'+row.protein+'.hmm', threads)
    os.remove(alignment)
    with open(hmm, 'r') as file:
        t1 = file.readlines()
        t1.insert(2, 'DESC  ' + row.desc + '\n')
    with open(dbpath + "irompdb/iromps.hmm", "a+") as f:
        f.writelines(t1)
    os.remove(hmm)
writer.save()
print('-' * int(get_terminal_size()[0]))
print("HMM library info written to:"+ dbpath +'irompdb/seed_alignment_blastp_' +
                        evalue + '_%i.xlsx' % window)
print("Done! Enjoy your new HMM library ;D")
