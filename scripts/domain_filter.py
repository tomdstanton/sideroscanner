#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'


def domain_filter(in_file, lowqual):
    q = run_hmmsearch(in_file, hmmpath+'/PF07715.hmm', threads)
    in_file_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    filter1 = ''
    for h in q.hit_keys:
        filter1 = filter1 + in_file_dict[h].format("fasta")
    print(f'Filtered {str(filter1.count(">"))} proteins with {q.description}')
    if lowqual is True or len(filter1) == 0:
        return filter1
    else:
        q = run_hmmsearch(filter1, hmmpath+'/PF00593.hmm', threads)
        filter2 = ''
        filter1_dict = to_dict(parse(StringIO(filter1), 'fasta'))
        for h in q.hit_keys:
            filter2 = filter2 + filter1_dict[h].format("fasta")
        print(f'Filtered {str(filter2.count(">"))} proteins with {q.description}')
        return filter2