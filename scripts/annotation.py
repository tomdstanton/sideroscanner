#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'


def annotation(in_file, input_type):
    data = []
    for q in run_hmmscan(in_file, parse_args().lib, threads):
        if len(q.hits) > 0:
            data.append(f'{q.id.rsplit("_", 1)[0]},{q.id}'
                        f',{q.hits[0].id},{q.hits[0].description}'
                        f',{str(q.hits[0].bitscore)}')
        else:
            return None
    print(f'Annotated {len(data)} proteins')
    hmmscan_df = pd.DataFrame([sub.split(",") for sub in data],
                              columns=['contig', 'query', 'hit', 'description', 'score'])
    data = []
    for r in parse(StringIO(in_file), 'fasta'):
        data.append(f'{r.id},{str(len(r.seq))},'
                    f'{str(round((ProteinAnalysis(r._seq._data).molecular_weight()/1000),1))}')
    len_mass_df = pd.DataFrame([sub.split(",") for sub in data],
                               columns=['query', 'len', 'kDa'])
    if not input_type == 'genome':
        return hmmscan_df.merge(len_mass_df, on='query', how='left')
    else:
        data = []
        for r in parse(StringIO(in_file), 'fasta'):
            data.append(','.join(r.description.replace(' ', '').split('#')[0:4]))
        bed_df = pd.DataFrame([sub.split(",") for sub in data],
                              columns=['query', 'start', 'end', 'str'])
        return hmmscan_df.merge(len_mass_df, on='query', how='left').merge(bed_df, on='query', how='left')