#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import zipfile
import gzip
from os import remove, uname
from argparse import RawTextHelpFormatter, ArgumentParser
from shutil import get_terminal_size
from datetime import datetime
from pathlib import Path
from Bio.SeqIO import parse, write

from scripts.fetch import fetch
from scripts.config import flankpath, mgepath, plspath
from scripts.blast import run_makeblastdb
from scripts.cdhit import run_cdhit

def parse_args():
    parser = ArgumentParser(add_help=False,
                                     formatter_class=RawTextHelpFormatter,
                                     usage="./builddbs.py -db [chosen_db]",
                                     description='''
        SideroScannner: a tool for annotating IROMPs in bacteria
        ========================================================
        Please cite: Stanton et al, 2020''')
    group = parser.add_argument_group("Options")
    group.add_argument('-db',  nargs='?', choices=['plsdb', 'mgedb', 'flankdb'],
                        help='''| choose a specific database to build
    [default: all]
-----------------------------------------------''')
    group.add_argument('-k', action='store_true',
                       help='''| keep database fasta files
    [default: all]
-----------------------------------------------''')

    group.add_argument("-h", action="help", help='''| show this help message and exit''')
    return parser.parse_args()

def get_plsdb():
    Path(plspath).mkdir(parents=True, exist_ok=True)
    fetch('https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip', None, plspath+'/plsdb')
    with zipfile.ZipFile(plspath+'/plsdb', "r") as zip_ref:
        zip_ref.extractall(plspath + '/')
    remove(plspath + '/plsdb')
    remove(plspath + '/plsdb.tsv')
    remove(plspath + '/plsdb_changes.tsv')
    remove(plspath + '/plsdb.msh')
    remove(plspath + '/plsdb.abr')
    remove(plspath + '/plsdb.sim')
    remove(plspath + '/README.md')
    return

def get_mgedb():
    Path(mgepath).mkdir(parents=True, exist_ok=True)
    ice = fetch('https://db-mml.sjtu.edu.cn/ICEberg2/download/ICE_seq_all.fas', None, mgepath + '/ice.fna')
    t4ss = fetch('https://db-mml.sjtu.edu.cn/ICEberg2/download/T4SS-type_ICE_seq_all.fas', None, mgepath + '/t4ss.fna')
    aice = fetch('https://db-mml.sjtu.edu.cn/ICEberg2/download/AICE_seq_all.fas', None, mgepath + '/aice.fna')
    ime = fetch('https://db-mml.sjtu.edu.cn/ICEberg2/download/IME_seq_all.fas', None, mgepath + '/ime.fna')
    cime = fetch('https://db-mml.sjtu.edu.cn/ICEberg2/download/CIME_seq_all.fas', None, mgepath + '/cime.fna')
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
    patric = fetch('ftp://ftp.patricbrc.org/specialty_genes/referenceDBs/PATRIC_VF.faa', None, flankpath + '/patric.faa')
    victors = fetch('http://www.phidias.us/victors/downloads/gen_downloads_protein.php', None, flankpath + '/victors.faa')
    vfdb = fetch('http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz', None, flankpath + '/vfdb.faa.gz')
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
    print('-' * int(get_terminal_size()[0]))
    print(__title__ + ': ' + __version__)
    print('Your system is ' + uname()[0])
    if 'Linux' not in uname()[0]:
        print('Warning: sideroscanner has not been tested on ' + uname()[0])
    print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))


    if parse_args().db is None:
        db_list = ['plsdb', 'mgedb', 'flankdb']
    else:
        db_list = parse_args().db

    if 'plsdb' in db_list:
        print('-' * int(get_terminal_size()[0]))
        print("Preparing plasmid database...")
        print('-' * int(get_terminal_size()[0]))
        get_plsdb()

    if 'mgedb' in db_list:
        print('-' * int(get_terminal_size()[0]))
        print("Preparing mobile genetic element database...")
        print('-' * int(get_terminal_size()[0]))
        mgedb = get_mgedb()
        print('-' * int(get_terminal_size()[0]))
        run_makeblastdb(mgedb, 'nucl', mgepath+'/mgedb')
        if parse_args().k is None:
            remove(mgedb)

    if 'flankdb' in db_list:
        print('-' * int(get_terminal_size()[0]))
        print("Preparing flanking virulence gene database...")
        print('-' * int(get_terminal_size()[0]))
        flankdb = get_flankdb()
        print('-' * int(get_terminal_size()[0]))
        flankdb = run_cdhit(flankdb, 1)
        print('-' * int(get_terminal_size()[0]))
        run_makeblastdb(flankdb, 'prot', flankpath+'/flankdb')
        if parse_args().k is None:
            remove(flankdb)

    print('-' * int(get_terminal_size()[0]))
    print("Done!")

if __name__ == "__main__":
    main()
