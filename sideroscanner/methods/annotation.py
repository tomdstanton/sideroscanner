#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sideroscanner.tools.hmmer3 import run_hmmscan
import pandas as pd
from Bio.SeqIO import to_dict, parse
from io import StringIO
from sideroscanner.tools.hmmer3 import run_hmmsearch


def domain_filter(in_file, iromppath, lowqual, threads):
    q = run_hmmsearch(in_file, iromppath + '/PF07715.hmm', threads)
    in_file_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    filter1 = ''
    for h in q.hit_keys:
        filter1 = filter1 + in_file_dict[h].format("fasta")
    print(f'Filtered {str(filter1.count(">"))} proteins with {q.description}')
    if lowqual is True or len(filter1) == 0:
        return filter1
    else:
        q = run_hmmsearch(filter1, iromppath + '/PF00593.hmm', threads)
        filter2 = ''
        filter1_dict = to_dict(parse(StringIO(filter1), 'fasta'))
        for h in q.hit_keys:
            filter2 = filter2 + filter1_dict[h].format("fasta")
        print(f'Filtered {str(filter2.count(">"))} proteins with {q.description}')
        return filter2


def annotate(in_file, input_type, lib, threads):
    data = []
    for q in run_hmmscan(in_file, lib, threads):
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
