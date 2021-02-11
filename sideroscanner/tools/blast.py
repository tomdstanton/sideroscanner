#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SearchIO
from subprocess import Popen, PIPE, DEVNULL, run
from io import StringIO
from sys import exit
from shutil import which
from os import path


def run_blastn(query, blastn_db, cov, percid, threads):
    if which('blastn') is None:
        exit('[!] blastn not found')
    print('[>] Running blastn...    ', end="", flush=True)
    cmd = ['blastn',
           '-db', blastn_db,
           '-task', 'megablast',
           # '-soft_masking', 'true',
           '-culling_limit', '1',
           '-outfmt', '6',
           # '-max_hsps', '1',
           # '-evalue', '1e-50',
           '-num_threads', threads,
           '-max_target_seqs', '1',
           '-perc_identity', str(percid),
           '-qcov_hsp_perc', str(cov),
           '-query', '-']
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    child.stdin.write(query.encode())
    return SearchIO.parse(StringIO(child.communicate()[0].decode('utf-8')), 'blast-tab')


def run_blastp(query, blastp_db, evalue, max_seqs, threads):
    if which('blastp') is None:
        exit('[!] blastp not found')
    print('[>] Running blastp...    ', end="", flush=True)
    cmd = ['blastp',
           '-db', blastp_db,
           '-outfmt', '5',
           #'-task', 'blastp-fast',
           '-evalue', evalue,
           # '-max_hsps', '1',
           '-max_target_seqs', max_seqs,
           # '-soft_masking', 'true',
           #'-word_size', '20',
           #'-threshold', '21',
           '-num_threads', threads,
           '-query', '-']
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    child.stdin.write(query.encode())
    return SearchIO.parse(StringIO(child.communicate()[0].decode('utf-8')), 'blast-xml')


def run_makeblastdb(infile, molecule, name):
    if which('makeblastdb') is None:
        exit('[!] makeblastdb not found')
    print('[>] Running makeblastdb...   ', end="", flush=True)
    cmd = ['makeblastdb', '-in', infile,
           '-title', name.rsplit('/', 1)[1].split('.')[0],
           '-parse_seqids',
           '-dbtype', molecule, '-out', name]
    run(cmd, stdout=DEVNULL)
    print('done!')

    return name
