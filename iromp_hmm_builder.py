#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import re
import pandas as pd
from os import cpu_count, remove, uname, path
from argparse import RawTextHelpFormatter, ArgumentParser
from Bio import SeqIO, Entrez
from shutil import get_terminal_size
from datetime import datetime

from scripts.blast import run_blastp, run_makeblastdb
from scripts.cdhit import run_cdhit
from scripts.fetch import fetch
from scripts.hmmer3 import run_hmmbuild
from scripts.trimal import run_trimal
from scripts.muscle import run_muscle
from scripts.config import iromppath, dbpath


def parse_args():
    parser = ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter,
                                     usage="./iromp_hmm_builder.py",
                                     description='''
            SideroScannner: a tool for annotating IROMPs in bacteria
            ========================================================''')
    group = parser.add_argument_group("Options")
    group.add_argument('-w', metavar='int', type=int, default=2,
                       help='''| protein length (w)indow for blastp
    [default: 3]
-----------------------------------------------''')
    group.add_argument('-e', metavar='str', type=str, default='1e-250',
                       help='''| (e)value for blastp
    [default: 1e-250]
-----------------------------------------------''')
    group.add_argument('-t', metavar='int', type=int, default=cpu_count(),
                       help='''| number of (t)hreads to use
    [default: %i]
-----------------------------------------------''' % cpu_count())
    group.add_argument("-h", action="help", help='''| show this help message and exit''')
    return parser.parse_args()


def get_fastadb():
    url = 'http://www.uniprot.org/uniprot/'
    params = {'query': 'taxonomy:"Bacteria [2]" AND '
                       'family:"tonb-dependent receptor family" AND '
                       'locations:(location:"Cell outer membrane [SL-0040]") '
                       'NOT partial NOT fragment', 'format': 'fasta'}
    print("Fetching TonB-depenent receptors")
    fetch(url, params, iromppath+'/iromps.faa')
    return iromppath+'/iromps.faa'


def seed_blast(in_file, fastadb, blastdb, name, db_df, excel):
    hit_acc = []
    for q in run_blastp(in_file, blastdb, parse_args().e, threads):
        if len(q.hits) > 0:
            print('%i hits found' % len(q.hits))
            up = q.seq_len + parse_args().w
            down = q.seq_len - parse_args().w
            for h in q.hits:
                if down <= h.seq_len <= up:
                    hit_acc.append(h.blast_id)
            print("%i hits filtered between %i-%i aa" % (len(hit_acc), down, up))
            hit_df = db_df[db_df['accession'].isin(hit_acc)]
            hit_df.to_excel(excel, sheet_name=name, index=False)
            hits = ''
            db_dict = SeqIO.to_dict(SeqIO.parse(fastadb, 'fasta'))
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

def main():
    print('-' * int(get_terminal_size()[0]))
    print(__title__ + ': ' + __version__)
    print('Your system is ' + uname()[0])
    if 'Linux' not in uname()[0]:
        print('Warning: sideroscanner has not been tested on ' + uname()[0])
    print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))
    print('Using ' + threads + ' threads...')
    if not path.isfile(iromppath + '/iromps_nr.faa'):
        if not path.isfile(iromppath + '/iromps.faa'):
            get_fastadb()
        else:
            fastadb = run_cdhit(iromppath + '/iromps.faa', 1)
    else:
        fastadb = iromppath+'/iromps_nr.faa'
    
    for i in [".psq", ".psi", ".psd", ".pog", ".pin", ".phr"]:
        if not path.isfile(iromppath+'/irompdb'+i):
            blastdb = run_makeblastdb(fastadb, 'prot', iromppath+'/irompdb')
            break
        else:
            blastdb = iromppath+'/irompdb'

    window = parse_args().w
    evalue = parse_args().e

    print("Window size of "+window+" and a blastp E-value of "+evalue)
    
    db_df = make_db_df(fastadb)
    seed_df = pd.read_csv(dbpath+'/iromps.csv')
    excel = pd.ExcelWriter(iromppath+'/seed_alignment_blastp_' +
                            evalue + '_%i.xlsx' % window, engine='xlsxwriter')
    # Iterate through each seed
    for row in seed_df.itertuples():
        # Download seed from NCBI
        print('-' * int(get_terminal_size()[0]))
        iromp = fetch_seed(row.acc, row.protein)

        # Blast against IROMP uniprot database
        hits = seed_blast(iromp, fastadb, blastdb, row.protein, db_df, excel)

        # Align hits with muscle
        print('-' * int(get_terminal_size()[0]))
        alignment = run_muscle(hits, iromppath + '/' + row.protein + '_aligned.faa')

        # Trim alignment with Trimal
        trimmed = run_trimal(alignment)

        # Make HMM from trimmed alignment
        print('-' * int(get_terminal_size()[0]))
        hmm = run_hmmbuild(trimmed, row.protein, iromppath+'/'+row.protein+'.hmm', threads)

        # Add descriptions to HMMs and concat them
        with open(hmm, 'r') as file:
            t1 = file.readlines()
            t1.insert(2, 'DESC  ' + row.desc + '\n')
        with open(iromppath+"/iromps.hmm", "a+") as f:
            f.writelines(t1)

        # Cleanup
        remove(hmm)
        remove(alignment)

    excel.save()
    print('-' * int(get_terminal_size()[0]))
    print("HMM library info written to:"+ iromppath +'/seed_alignment_blastp_' +
                            evalue + '_%i.xlsx' % window)
    print("Done! Enjoy your new HMM library ;D")
    
if __name__ == "__main__":
    if parse_args().t > cpu_count():
        print('Number of threads exceeds available CPUs, will use: %i' % cpu_count())
    else:
        threads = str(parse_args().t)
    main()
