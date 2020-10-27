#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from zipfile import ZipFile
from os import remove

from sideroscanner.methods.fetch import fetch_url

def get_plsdb(plspath):
    print("[-] Preparing plasmid database")
    Path(plspath).mkdir(parents=True, exist_ok=True)
    fetch_url('https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip', None, plspath + '/plsdb')
    with ZipFile(f'{plspath}/plsdb', "r") as zip_ref:
        zip_ref.extractall(f'{plspath}/')
    for f in ['/plsdb.tsv','/plsdb_changes.tsv','/plsdb.msh','/plsdb.abr','/plsdb.sim','/README.md']:
        remove(plspath + f)
    return