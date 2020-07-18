#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import RawTextHelpFormatter, ArgumentParser
from datetime import datetime
from os import cpu_count, remove, uname, path
from os.path import dirname, abspath
from shutil import get_terminal_size
from sys import argv

import pandas as pd
from Bio.SeqIO import to_dict, parse

from scripts.efetch import efetch
from scripts.get_fastadb import get_fastadb
from scripts.make_db_df import make_db_df
from scripts.seed_blast import seed_blast
from scripts.tools.blast import run_makeblastdb
from scripts.tools.cdhit import run_cdhit
from scripts.tools.hmmer3 import run_hmmbuild
from scripts.tools.muscle import run_muscle
from scripts.tools.trimal import run_trimal

dbpath = abspath(dirname(argv[0]))+'/data'
iromppath = dbpath+'/irompdb'
blastdb = iromppath+'/irompdb'
lib = iromppath+"/iromps.hmm"
proteins = iromppath+"/iromps.faa"
clustered = iromppath+"/iromps_nr.faa"


def parse_args():
    parser = ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter,
                                     usage="./sideroscanner-buildhmms.py",
                                     description='''
SideroScannner: a tool for annotating IROMPs in bacteria
----------------------------------------------------------
Please cite: Stanton et al, 2020''')
    group = parser.add_argument_group("Options")
    group.add_argument('-o', metavar='-', nargs='?', type=str,
                       const=f'seed_alignment_{datetime.today().strftime("%d%m%y_%H%M")}.xlsx',
                       help='''save info about proteins in seed alignment
[optional: path/to/output/file]
[default: seed_alignment_DDMMYY_hhmm.xlsx]
-----------------------------------------------''')
    group.add_argument('-w', metavar='int', type=int, default=2,
                       help='''protein length window [default: 2]
-----------------------------------------------''')
    group.add_argument('-e', metavar='str', type=str, default='1e-130',
                       help='''evalue for blastp [default: 1e-130]
-----------------------------------------------''')
    group.add_argument('-t', metavar='int', type=int, default=cpu_count(),
                       help=f'''number of threads [default: {cpu_count()}]
-----------------------------------------------''')
    group.add_argument('--keep', action='store_true',default=False,
                       help='''keep seed alignments
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''show this help message and exit''')
    return parser.parse_args()


def main():
    print('-' * int(get_terminal_size()[0]))
    print(f'{__name__} {__version__}')
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

    append = False
    if path.isfile(lib):
        answer = None
        while answer not in ("append", "overwrite", "exit" ):
            answer = input(f'[!] IROMP HMM library already exists, append, overwrite or exit? [a/o/e]: ')
            if answer == "a":
                append = True
                break
            elif answer == "o":
                break
            elif answer == "e":
                exit('[!] Exiting')
            else:
                print("[!] Please enter either [a], [o] or [e]: ")

    if append is False:
        for i in [".h3i", ".h3p", ".h3f", ".h3m", ""]:
            if path.isfile(f'{lib}{i}'):
                remove(f'{lib}{i}')
    else:
        for i in [".h3i", ".h3p", ".h3f", ".h3m"]:
            if path.isfile(f'{lib}{i}'):
                remove(f'{lib}{i}')

    seed_df = pd.read_csv(iromppath+'/iromps.csv')
    total_seeds = len(seed_df.index)
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
    writer = pd.ExcelWriter(f'{parse_args().o}', engine='xlsxwriter')
    c = 0
    d = {}
    for row in seed_df.itertuples():
        if append is True:
            with open(lib) as f:
                if row.acc in f.read():
                    continue
        print('-' * int(get_terminal_size()[0]))
        c = c + 1
        print(f'[-] Building HMM: \033[4m{c} of {total_seeds}\033[0m')
        hit_acc = seed_blast(efetch(row.acc), blastdb, parse_args().e, parse_args().w, threads)
        if hit_acc is None:
            print(f'[!] No hits found for {row.protein} ({row.acc}), review entry in iromps.csv')
            continue
        hit_df = db_df[db_df['accession'].isin(hit_acc)]
        hit_df.to_excel(writer, sheet_name=row.protein, index=False)
        hits = ''
        for h in hit_acc:
            hits = hits + db_dict[h].format("fasta")
        d[row.protein] = str(hit_acc)
        alignment = run_muscle(hits, f'{iromppath}/{row.protein}_aligned.faa')
        trimmed = run_trimal(alignment)
        if parse_args().keep is False:
            remove(alignment)
        hmm = run_hmmbuild(trimmed, row.protein, f'{iromppath}/{row.protein}.hmm', threads)
        with open(hmm, 'r') as file, open(lib, "a+") as f:
            t1 = file.readlines()
            t1.insert(2, 'ACC   ' + row.acc + '\n')
            t1.insert(3, 'DESC  ' + row.desc + '\n')
            f.writelines(t1)
        print(f'added {row.protein} to HMM library')
        remove(hmm)
    print('-' * int(get_terminal_size()[0]))
    rev_d = {}
    for key, value in d.items():
        try:
            rev_d[value].append(key)
        except:
            rev_d[value] = [key]
    overlap = [value for key, value in rev_d.items() if len(value) > 1]
    if len(overlap) > 0:
        print(f'[!] {", ".join(overlap[0])} share proteins in their seed alignments, revise entries')
        print('-' * int(get_terminal_size()[0]))
    print(f'[-] HMMs written to: {iromppath}/iromps.hmm')
    if parse_args().o is not None:
        writer.save()
        print(f'[-] Info written to: {parse_args().o}')

if __name__ == "__main__":
    if parse_args().t > cpu_count():
        print('[!] Number of threads exceeds available CPUs, will use: %i' % cpu_count())
    else:
        threads = str(parse_args().t)
    main()
