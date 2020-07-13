#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'


def export_proteins(in_file, hits, seq):
    to_write = []
    in_file_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    for row in hits.itertuples():
        r = in_file_dict[row.query]
        r.id = row.query
        r.description = f'{row.description} {row.hit}'
        to_write.append(r)
    if parse_args().e == 'seq':
        filename = f'{path.splitext(seq)[0]}_sideroscanner.faa'
    else:
        filename = parse_args().e
    write(to_write, filename, "fasta")
    print(f'Proteins written to: {filename}')