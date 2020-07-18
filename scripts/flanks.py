#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.SeqIO import parse, to_dict
from io import StringIO
import re
from .tools.blast import run_blastp
import pandas as pd


def flank_screen(in_file, hits, flankpath, cds, threads):
    cds_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    queries = ''
    for h in hits['query'].tolist():
        pos = int(h.split('_')[-1])
        acc = f'{h.rsplit("_", 1)[0]}_'
        for i in range(cds):
            i = i + 1
            try:
                u = cds_dict[acc + str(pos - i)]
                u.id = h
                u.description = f'up_{str(i)}'
                queries = queries + u.format("fasta")
            except:
                continue
            try:
                d = cds_dict[acc + str(pos + i)]
                d.id = h
                d.description = f'down_{str(i)}'
                queries = queries + d.format("fasta")
            except:
                continue
    up = []
    down = []
    for q in run_blastp(queries, flankpath+'/flankdb', '1e-130', threads):
        if len(q.hits) > 0:
            gene_name = re.sub("\s*\[[^[]*\]$", '', q.hsps[0].hit_description).strip()
            if gene_name.startswith('('):
                gene_name = gene_name.split('(', 1)[1].split(')')[0]
            if 'up' in q.hsps[0].query_description:
                up.append(f'{q.hsps[0].query_id}#{gene_name}')
            if 'down' in q.hsps[0].query_description:
                down.append(f'{q.hsps[0].query_id}#{gene_name}')
    print(f'{len(up)} upstream CDS and {len(down)} downstream CDS identified')
    u_df = pd.DataFrame([sub.split("#") for sub in up], columns=['query', 'upstream'])
    d_df = pd.DataFrame([sub.split("#") for sub in down], columns=['query', 'downstream'])
    u_df = u_df.groupby(['query'])['upstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    d_df = d_df.groupby(['query'])['downstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    return hits.merge(u_df, on='query', how='left').merge(d_df, on='query', how='left')