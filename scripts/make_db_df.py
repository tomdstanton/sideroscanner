#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.SeqIO import parse
from re import search
import pandas as pd

def make_db_df(in_file):
    print('[>] Creating DataFrame... ', end="", flush=True)
    data = []
    for r in parse(in_file, 'fasta'):
        if search(r'GN=(.*) PE=', r.description):
            gene = search(r'GN=(.*) PE=', r.description).group(1)
        else:
            gene = 'NA'
        data.append(f'{r.id}${gene}$'
                    f'{search(r" (.*) OS=", r.description).group(1).replace(",", " ")}$'
                    f'{search(r"OS=(.*) OX=", r.description).group(1).replace(",", " ")}')
    print(f'{len(data)} total proteins in database')
    return pd.DataFrame([sub.split("$") for sub in data],
                        columns=['accession', 'gene', 'description', 'organism'])