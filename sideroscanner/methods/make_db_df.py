#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from re import search
import pandas as pd

def make_db_df(in_file):
    print('[>] Creating DataFrame...    ', end="", flush=True)

    with open(in_file) as f:
        headers = [line.rstrip().replace('>', '') for line in f if line.startswith('>')]

    data = []
    for header in headers:
        accession = header.split(' ', 1)[0]
        description = header.split(' ', 1)[1]

        if search(r'GN=(.*) PE=', description):
            gene = search(r'GN=(.*) PE=', description).group(1)
        else:
            gene = 'NA'

        info = description.split(' OS=')[0]
        org = description.split(' OS=')[1].split(' OX=')[0]

        data.append(f'{accession}${gene}${org}${info}')

    print(f'{len(data)} total proteins in database')
    return pd.DataFrame([sub.split("$") for sub in data],
                        columns=['accession', 'gene', 'description', 'organism'])