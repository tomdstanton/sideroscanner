#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from argparse import RawTextHelpFormatter, ArgumentParser
from datetime import datetime
from os import cpu_count, remove, uname, path
from os.path import dirname, abspath
from shutil import get_terminal_size
from sys import argv

import pandas as pd
from Bio.SeqIO import to_dict, parse

from .scripts.efetch import efetch
from .scripts.get_fastadb import get_fastadb
from .scripts.make_db_df import make_db_df
from .scripts.seed_blast import seed_blast
from .scripts.tools.blast import run_makeblastdb
from .scripts.tools.cdhit import run_cdhit
from .scripts.tools.hmmer3 import run_hmmbuild
from .scripts.tools.muscle import run_muscle
from .scripts.tools.trimal import run_trimal

dbpath = abspath(dirname(argv[0]))+'/data'
iromppath = dbpath+'/irompdb'
blastdb = iromppath+'/irompdb'
lib = iromppath+"/iromps.hmm"
proteins = iromppath+"/iromps.faa"
clustered = iromppath+"/iromps_nr.faa"


def parse_args():
    parser = ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter,
                                     usage="./scripts-buildhmms.py",
                                     description='''
            SideroScannner: a tool for annotating IROMPs in bacteria
            ========================================================''')
    group = parser.add_argument_group("Options")
    group.add_argument('-w', metavar='int', type=int, default=2,
                       help='''protein length window [default: 2]
-----------------------------------------------''')
    group.add_argument('-e', metavar='str', type=str, default='1e-130',
                       help='''evalue for blastp [default: 1e-130]
-----------------------------------------------''')
    group.add_argument('-t', metavar='int', type=int, default=cpu_count(),
                       help=f'''number of threads [default: {cpu_count()}]
-----------------------------------------------''')
    group.add_argument('--keep', action='store_true',
                       help='''keep seed alignments and IROMP database files
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''show this help message and exit''')
    return parser.parse_args()


def main():
    print('-' * int(get_terminal_size()[0]))
    print(f'{__title__} {__version__}')
    print(f'Your system is {uname()[0]}')
    if 'Linux' not in uname()[0]:
        print(f'[!] Warning: scripts has not been tested on {uname()[0]}')
    print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))
    print(f'Using {threads} threads...')
    print('-' * int(get_terminal_size()[0]))

    if not path.isdir(iromppath):
        exit(f'[!] {iromppath} is not an existing directory')

    if not path.isfile(iromppath+'/iromps.csv'):
        exit(f'[!] {iromppath}/iromps.csv not found')

    seed_df = pd.read_csv(iromppath+'/iromps.csv')
    total_seeds = len(seed_df.index)

    answer = None
    while answer not in ("y", "n"):
        answer = input(f'[!] Build pHMMs for {total_seeds} IROMPs '
          f'using an E-value of {parse_args().e} and '
          f'length window of {parse_args().w} [y/n]: ')
        if answer == "y":
            continue
        elif answer == "n":
            exit('[!] Exiting')
        else:
            print("[!] Please enter y or n: ")

    files = [] 
    for f in [proteins, clustered, blastdb+'.psq']:
        files.append(path.isfile(f))
    
    if True not in files:
        fastadb = run_cdhit(get_fastadb(iromppath), 1)
        run_makeblastdb(fastadb, 'prot', blastdb)

    elif files[1] is False and files[2] is False:
        fastadb = run_cdhit(proteins, 1)
        run_makeblastdb(fastadb, 'prot', blastdb)

    elif files[2] is False:
        run_makeblastdb(clustered, 'prot', blastdb)
        fastadb = clustered
    else:
        fastadb = clustered

    db_df = make_db_df(fastadb)
    db_dict = to_dict(parse(fastadb, 'fasta'))
    writer = pd.ExcelWriter(f'{iromppath}/seed_alignment_blastp_{parse_args().e}_{parse_args().w}.xlsx',
                                engine='xlsxwriter')

    counter = 0
    # Iterate through each seed
    for row in seed_df.itertuples():
        print('-' * int(get_terminal_size()[0]))
        counter = counter + 1
        print(f'[-] Building pHMM: \033[4m{counter} of {total_seeds}\033[0m')
        hit_acc = seed_blast(efetch(row.acc), blastdb, parse_args().e, parse_args().w, threads)
        hit_df = db_df[db_df['accession'].isin(hit_acc)]
        hit_df.to_excel(writer, sheet_name=row.protein, index=False)
        hits = ''
        for h in hit_acc:
            hits = hits + db_dict[h].format("fasta")
        alignment = run_muscle(hits, f'{iromppath}/{row.protein}_aligned.faa')
        trimmed = run_trimal(alignment)
        remove(alignment)
        hmm = run_hmmbuild(trimmed, row.protein, f'{iromppath}/{row.protein}.hmm', threads)
        with open(hmm, 'r') as file:
            t1 = file.readlines()
            t1.insert(2, 'DESC  ' + row.desc + '\n')
        with open(lib, "a+") as f:
            f.writelines(t1)
        print(f'added {row.protein} to HMM library')
        remove(hmm)

    writer.save()
    print('-' * int(get_terminal_size()[0]))
    print(f'[-] HMMs written to: {iromppath}/iromps.hmm')
    print(f'[-] Info written to: {iromppath}/seed_alignment_blastp_'
          f'{parse_args().e}_{parse_args().w}.xlsx')

if __name__ == "__main__":
    if parse_args().t > cpu_count():
        print('[!] Number of threads exceeds available CPUs, will use: %i' % cpu_count())
    else:
        threads = str(parse_args().t)
    main()
