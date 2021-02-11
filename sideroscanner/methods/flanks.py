#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.SeqIO import parse, to_dict
from io import StringIO
from sideroscanner.tools.blast import run_blastp
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
                if len(u.seq) > 2000:
                    print(f'[!] Warning: {u.id} is >2000aa and might cause blastp to hang')
                u.id = h
                u.description = f'{i} upstream'
                queries = queries + u.format("fasta")
            except:
                continue
            try:
                d = cds_dict[acc + str(pos + i)]
                if len(d.seq) > 2000:
                    print(f'[!] Warning: {d.id} is >2000aa and might cause blastp to hang')
                d.id = h
                d.description = f'{i} downstream'
                queries = queries + d.format("fasta")
            except:
                continue

    hit_list = []
    # For each IROMP flanking protein
    for q in run_blastp(queries, flankpath+'/flankdb', '1e-50', '5', threads):
        # If the protein has a hit
        if len(q.hits) > 0:
            out = q.hsps[0].hit_description
            if out.startswith('('):
                gene_name = out.split('(', 1)[1].split(')')[0]
                description = out.split(')', 1)[1].split('[', 1)[0]
            elif '[' in out:
                gene_name = ''
                description = out.split('[', 1)[0]
            else:
                gene_name = out.split('GN=', 1)[1].split('PE=')[0]
                description = out.split('OS=')[0]
            hit_list.append(f'{q.hsps[0].query_id}#'
                            f'{q.hsps[0].query_description}: '
                            f'{gene_name.strip()} '
                            f'{description.strip()}\n')

    print(f'{len(hit_list)} flanking genes identified')

    hit_df = pd.DataFrame([sub.split("#") for sub in hit_list],
                        columns=['query','flanking_genes'])

    hit_df = hit_df.groupby(['query'])['flanking_genes'].apply(lambda x: ''.join(x.astype(str))).reset_index()
    hit_df['flanking_genes'] = hit_df['flanking_genes'].str.strip()
    return hits.merge(hit_df, on='query', how='left')