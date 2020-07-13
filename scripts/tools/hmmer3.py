#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from Bio import SearchIO
from subprocess import Popen, PIPE, run
from io import StringIO

def run_hmmsearch(in_file, hmm, threads):
    print("[>] Running hmmsearch... ", end="", flush=True)
    cmd = ['hmmsearch', '--noali', '--cut_tc', '--cpu', threads, hmm, '-']
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    child.stdin.write(in_file.encode())
    return SearchIO.read(StringIO(child.communicate()[0].decode('utf-8')), 'hmmer3-text')


def run_hmmscan(in_file, hmm, threads):
    print("[>] Running hmmscan...   ", end="", flush=True)
    cmd = ['hmmscan', '--noali',
           #'--domE', '0.2', '-E', '1e-50',
           '--cpu', threads, hmm, '-']
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    child.stdin.write(in_file.encode())
    return SearchIO.parse(StringIO(child.communicate()[0].decode('utf-8')), 'hmmer3-text')


def run_hmmbuild(in_file, id, out_file, threads):
    print("[>] Running hmmbuild...")
    cmd = ['hmmbuild', '--cpu', threads,
           '--informat', 'clustal', '--amino',
           '-n', id, out_file, '-']
    child = Popen(cmd, stdin=PIPE)
    child.stdin.write(in_file.encode())
    child.communicate()
    return out_file


def run_hmmpress(lib):
    run(["hmmpress", lib])