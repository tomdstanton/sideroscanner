#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import zipfile
from argparse import RawTextHelpFormatter, ArgumentParser
from datetime import datetime
from os import remove, uname
from os.path import dirname, abspath, isfile
from glob import glob
from pathlib import Path
from shutil import get_terminal_size
from sys import argv

from Bio.SeqIO import parse, write

from scripts.fetch import fetch_url
from scripts.tools.blast import run_makeblastdb
from scripts.tools.cdhit import run_cdhit

dbpath = abspath(dirname(argv[0]))+'/data'
mgepath = dbpath+'/mgedb'
plspath = dbpath+'/plsdb'
flankpath = dbpath+'/flankdb'

def parse_args():
    parser = ArgumentParser(add_help=False,
                                     formatter_class=RawTextHelpFormatter,
                                     usage="./sideroscanner-builddbs.py",
                                     description='''
SideroScannner: a tool for annotating IROMPs in bacteria
-----------------------------------------------------------------------
Please cite: Stanton et al, 2020''')
    group = parser.add_argument_group("Options")
    group.add_argument('-db', nargs='+',
                       choices=['plsdb', 'mgedb', 'flankdb'],
                       default=['plsdb', 'mgedb', 'flankdb'],
                       metavar='plsdb/mgedb/flankdb',
                        help='''build specific database [default: all]
-----------------------------------------------''')
    group.add_argument('--keep', action='store_true',
                       help='''keep database fasta files
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''show this help message and exit''')
    return parser.parse_args()


def get_plsdb():
    Path(plspath).mkdir(parents=True, exist_ok=True)
    fetch_url('https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip', None, plspath + '/plsdb')
    with zipfile.ZipFile(plspath+'/plsdb', "r") as zip_ref:
        zip_ref.extractall(plspath + '/')
    for f in ['/plsdb','/plsdb.tsv','/plsdb_changes.tsv','/plsdb.msh','/plsdb.abr','/plsdb.sim','/README.md']:
        remove(plspath + f)
    return

def get_mgedb():
    Path(mgepath).mkdir(parents=True, exist_ok=True)
    ice = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/ICE_seq_all.fas', None, mgepath + '/ice.fna')
    t4ss = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/T4SS-type_ICE_seq_all.fas', None, mgepath + '/t4ss.fna')
    aice = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/AICE_seq_all.fas', None, mgepath + '/aice.fna')
    ime = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/IME_seq_all.fas', None, mgepath + '/ime.fna')
    cime = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/CIME_seq_all.fas', None, mgepath + '/cime.fna')
    filenames = [ice, t4ss, aice, ime, cime]
    with open(mgepath + '/mgedb', 'wb') as mge_file:
        for fname in filenames:
            with open(fname, 'rb') as infile:
                for line in infile:
                    mge_file.write(line)
            remove(fname)
    with open(mgepath + '/mgedb_clean', 'a') as f:
        record_ids = list()
        for r in parse(mgepath + '/mgedb', 'fasta'):
            if len(r.id) > 50:
                r.id = r.id[0:49]
            if r.id not in record_ids:
                record_ids.append(r.id)
                write(r, f, 'fasta')
    remove(mgepath + '/mgedb')
    return mgepath + '/mgedb_clean'

def get_flankdb():
    Path(plspath).mkdir(parents=True, exist_ok=True)
    patric = fetch_url('ftp://ftp.patricbrc.org/specialty_genes/referenceDBs/PATRIC_VF.faa', None, flankpath + '/patric.faa')
    victors = fetch_url('http://www.phidias.us/victors/downloads/gen_downloads_protein.php', None, flankpath + '/victors.faa')
    vfdb = fetch_url('http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz', None, flankpath + '/vfdb.faa.gz')
    filenames = [patric, victors, vfdb]
    with open(flankpath + '/flankdb', 'wb') as flank_file:
        for fname in filenames:
            if fname.endswith('.gz'):
                with gzip.open(fname) as infile:
                    for line in infile:
                        flank_file.write(line)
            else:
                with open(fname, 'rb') as infile:
                    for line in infile:
                        flank_file.write(line)
            remove(fname)
    return flankpath + '/flankdb'

def main():
    if len(parse_args().db) == 0:
        exit('[!] -db flag requires either [plsdb/mgedb/flankdb]')

    print('-' * int(get_terminal_size()[0]))
    print(f'{__name__} {__version__}')
    print(f'Your system is {uname()[0]}')
    if 'Linux' not in uname()[0]:
        print(f'[!] Warning: scripts has not been tested on {uname()[0]}')
    print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))

    db_list = []
    for i in parse_args().db:
        if [n for n in glob(f'{dbpath}/{i}/{i}.*sq') if isfile(n)]:
            answer = None
            while answer not in ("y", "n"):
                answer = input(f'[!] {i} already exists, continue? [y/n]: ')
                if answer == "y":
                    db_list.append(i)
                elif answer == "n":
                    continue
                else:
                    print("[!] Please enter [y] or [n]: ")
        else:
            db_list.append(i)

    if not db_list:
        print('-' * int(get_terminal_size()[0]))
        exit('[!] All database files exist, exiting')

    if 'plsdb' in db_list:
        print('-' * int(get_terminal_size()[0]))
        print("[-] \033[4mPreparing plasmid database\033[0m")
        get_plsdb()

    if 'mgedb' in db_list:
        print('-' * int(get_terminal_size()[0]))
        print("[-] \033[4mPreparing mobile genetic element database\033[0m")
        mgedb = get_mgedb()
        if parse_args().keep is None:
            remove(mgedb)

    if 'flankdb' in db_list:
        print('-' * int(get_terminal_size()[0]))
        print("[-] \033[4mPreparing flanking virulence gene database\033[0m")
        flankdb = get_flankdb()
        flankdb = run_cdhit(flankdb, 1)
        run_makeblastdb(flankdb, 'prot', flankpath+'/flankdb')
        if parse_args().keep is None:
            remove(flankdb)

if __name__ == "__main__":
    main()
