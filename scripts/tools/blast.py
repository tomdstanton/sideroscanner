#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from Bio import SearchIO
from subprocess import Popen, PIPE, run
from io import StringIO

def run_blastn(query, blastn_db, cov, percid, threads):
    print('[>] Running blastn...    ', end="", flush=True)
    cmd = ['blastn',
           '-query', query,
           '-db', blastn_db,
           '-soft_masking', 'true',
           '-culling_limit', '1',
           '-outfmt', '6',
           '-max_hsps', '1',
           '-evalue', '1e-50',
           '-num_threads', threads,
           '-max_target_seqs', '1',
           '-perc_identity', str(percid),
           '-qcov_hsp_perc', str(cov)]
    child = Popen(cmd, stdout=PIPE, stderr=PIPE)
    return SearchIO.parse(StringIO(child.communicate()[0].decode('utf-8')), 'blast-tab')


def run_blastp(query, blastp_db, evalue, threads):
    print('[>] Running blastp...    ', end="", flush=True)
    cmd = ['blastp',
           '-db', blastp_db,
           '-outfmt', '5',
           '-max_hsps', '1',
           '-evalue', evalue,
           '-max_target_seqs', '100000',
           '-soft_masking', 'true',
           '-num_threads', threads,
           '-query', '-']
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    child.stdin.write(query.encode())
    return SearchIO.parse(StringIO(child.communicate()[0].decode('utf-8')), 'blast-xml')

def run_makeblastdb(infile, molecule, name):
    print('[>] Running makeblastdb...   ', end="", flush=True)
    cmd = ['makeblastdb', '-in', infile,
           '-title', name.rsplit('/', 1)[1].split('.')[0],
           '-parse_seqids',
           '-dbtype', molecule, '-out', name]
    run(cmd)
    return name