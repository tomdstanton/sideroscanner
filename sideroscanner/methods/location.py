#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sideroscanner.tools.blast import run_blastn

def range_subset(range1, range2):
    if not range1:
        return True
    if not range2:
        return False
    if len(range1) > 1 and range1.step % range2.step:
        return False
    return range1.start in range2 and range1[-1] in range2


def plasmid_mge_screen(in_file, plspath, mgepath, hits, percid, threads):
    plasmids = []
    for q in run_blastn(in_file, plspath+'/plsdb.fna', 10, percid, threads):
        for h in q.hits:
            plasmids.append(f'{h.query_id}#{h.id}')
    print(f'{len(plasmids)} plasmid(s) found')
    mges = []
    for q in run_blastn(in_file, mgepath+'/mgedb', 0, percid, threads):
        for h in q.hits:
            mges.append(f'{h.query_id}#{h.id.split("|")[2]}#'
                        f'{str(h.hsps[0].query_range[0])}-'
                        f'{str(h.hsps[0].query_range[1])}')
    print("%i MGE(s) found"%len(mges))
    for row in hits.itertuples():
        for line in mges:
            start = int(line.split('#')[2].split('-')[0])
            end = int(line.split('#')[2].split('-')[1])
            if range_subset(range(int(row.start), int(row.end)), range(start, end)):
                hits.at[row.Index, 'plasmid/mge'] = line.split('#')[1]
    for row in hits.itertuples():
        for line in plasmids:
            if row.contig == line.split('#')[0]:
                hits.at[row.Index, 'plasmid/mge'] = line.split('#')[1]
    return hits