#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from io import StringIO
from Bio.SeqIO import write, parse, to_dict

def export_proteins(in_file, out_file, name, hits):
    to_write = []
    in_file_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    for row in hits.itertuples():
        r = in_file_dict[row.query]
        r.id = row.query
        r.description = f'{row.description} {row.hit}'
        to_write.append(r)
    if out_file == '':
        if name == '-':
            out_file = 'stdin_sideroscanner.faa'
        else:
            out_file = f'{name}_sideroscanner.faa'
    write(to_write, out_file, "fasta")
    print(f'[-] Proteins written to: {out_file}')