#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from io import StringIO
from Bio.SeqIO import write, parse, to_dict
from os.path import splitext

def export_proteins(in_file, filename, hits, seq):
    to_write = []
    in_file_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    for row in hits.itertuples():
        r = in_file_dict[row.query]
        r.id = row.query
        r.description = f'{row.description} {row.hit}'
        to_write.append(r)
    if filename == 'seq':
        filename = f'{splitext(seq)[0]}_sideroscanner.faa'
    write(to_write, filename, "fasta")
    print(f'[-] Proteins written to: {filename}')
