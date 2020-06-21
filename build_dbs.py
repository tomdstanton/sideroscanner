#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import zipfile
import sys
import os
import argparse
from argparse import RawTextHelpFormatter

from scripts.fetch import fetch

pathname = os.path.dirname(sys.argv[0])
full_path = os.path.abspath(pathname)

def parse_args():
    parser = argparse.ArgumentParser(add_help=False,
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

dbpath = parse_args().dbpath


def get_mgedb():
    mge_file = dbpath + 'mgedb/mgedb'
    with open(mge_file, mode='wb') as localfile:
        localfile.write(fetch('http://202.120.12.136/ICEberg2/download/ICE_aa_all.fas'), None)


def get_plsdb():
    plsdb_file = dbpath + 'plsdb/plsdb'
    with open(plsdb_file, mode='wb') as localfile:
        localfile.write(fetch('https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip'), None)
    with zipfile.ZipFile(dbpath + 'plsdb', "r") as zip_ref:
        zip_ref.extractall(dbpath + 'plsdb/')
    os.remove(dbpath + 'plsdb')
    os.remove(dbpath + 'plsdb.tsv')
    os.remove(dbpath + 'plsdb_changes.tsv')
    os.remove(dbpath + 'plsdb.msh')
    os.remove(dbpath + 'plsdb.abr')
    os.remove(dbpath + 'plsdb.sim')
    os.remove(dbpath + 'README.md')

