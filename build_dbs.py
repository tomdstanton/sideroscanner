#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import zipfile
from os import cpu_count, remove
from argparse import RawTextHelpFormatter, ArgumentParser

from scripts.fetch import fetch
from scripts.config import full_path

def parse_args():
    parser = ArgumentParser(add_help=False,
                                     formatter_class=RawTextHelpFormatter,
                                     usage="./build_dbs.py --dbpath [path/]",
                                     description='''
        SideroScannner: a tool for annotating IROMPs in bacteria
        ========================================================''')
    group = parser.add_argument_group("Options")
    group.add_argument('--dbpath', metavar='path', type=str, default=full_path + '/databases/',
                        help='''| path/to/db/
    [default: ''' + full_path + '''/databases/]
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''| show this help message and exit''')
    return parser.parse_args()


def get_mgedb():
    mge_file = dbpath + '/mgedb'
    with open(mge_file, mode='wb') as file:
        file.write(fetch('http://202.120.12.136/ICEberg2/download/ICE_aa_all.fas'), None)


def get_plsdb():
    plsdb_file = dbpath + 'plsdb/plsdb'
    with open(plsdb_file, mode='wb') as localfile:
        localfile.write(fetch('https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip'), None)
    with zipfile.ZipFile(dbpath + 'plsdb', "r") as zip_ref:
        zip_ref.extractall(dbpath + 'plsdb/')
    remove(dbpath + 'plsdb')
    remove(dbpath + 'plsdb.tsv')
    remove(dbpath + 'plsdb_changes.tsv')
    remove(dbpath + 'plsdb.msh')
    remove(dbpath + 'plsdb.abr')
    remove(dbpath + 'plsdb.sim')
    remove(dbpath + 'README.md')

def main():
    xhd

if __name__ == "__main__":
    dbpath = parse_args().dbpath
    if parse_args().t > cpu_count():
        print('Number of threads exceeds available CPUs, will use: %i' % cpu_count())
    else:
        threads = str(parse_args().t)
    main()







