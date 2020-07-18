#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.SeqIO import parse, to_dict
import pandas as pd
from .tools.meme import run_mast

def tfbs_screen(in_file, furpath, hits, length):
    prom_dict = to_dict(parse(in_file, 'fasta'))
    queries = ''
    for row in hits.itertuples():
        if row.str == '1':
            f_start = int(row.start) - length
            f_end = int(row.start)
        elif row.str == '-1':
            f_start = int(row.end)
            f_end = int(row.end) + length
        try:
            x = prom_dict[row.contig][f_start:f_end]
            x.id = f'{row.query}:{str(f_start)}'
            queries = queries + x.format("fasta")
        except:
            continue

    mast_out = run_mast(queries, furpath+'/fur.meme').rstrip()
    results = []
    for line in (line for line in mast_out.split('\n') if not line.startswith('#')):
        results.append(line)
    if len(results) == 0:
        print("No binding sites found")
        return hits
    else:
        print(f'Putative TFBS for {len(results)} hit(s)')
        mast = []
        for q in results:
            query = q.split(' ')[0].split(':')[0]
            fur_start = int(q.split(' ')[0].split(':')[1]) + int(q.split(' ')[4])
            fur_end = int(q.split(' ')[0].split(':')[1]) + int(q.split(' ')[5])
            pval = str(q.split(' ')[8])
            bs = prom_dict[query.rsplit('_',1)[0]][fur_start:fur_end].seq._data
            if q.split(' ')[1] == '-1':
                bs = str(Seq(bs).reverse_complement())
            mast.append(f'{query}#{str(fur_start)}#{str(fur_end)}#{pval}#{bs}')
        df = pd.DataFrame([sub.split("#") for sub in mast],
                          columns=['query','fur_start','fur_end',
                                   'p_value','fur_box'])
        return hits.merge(df, on='query', how='left')