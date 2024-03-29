#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SearchIO
from subprocess import Popen, PIPE, DEVNULL
from io import StringIO
from shutil import which
from sys import exit

def run_hmmsearch(in_file, hmm, threads):
    if which('hmmsearch') is None:
        exit('[!] hmmsearch not found')
    print("[>] Running hmmsearch... ", end="", flush=True)
    cmd = ['hmmsearch', '--noali', '--cut_tc', '--cpu', threads, hmm, '-']
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    try:
        child.stdin.write(in_file.encode())
        return SearchIO.read(StringIO(child.communicate()[0].decode('utf-8')), 'hmmer3-text')
    except IOError:
        return print(f'\n[!] hmmer failed to run, is your input DNA?')


def run_hmmscan(in_file, hmm, threads):
    if which('hmmscan') is None:
        exit('[!] hmmscan not found')
    cmd = ['hmmscan', '--noali',
           #'--domE', '0.2', '-E', '1e-50',
           '--cpu', threads, hmm, '-']
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    try:
        child.stdin.write(in_file.encode())
        return SearchIO.parse(StringIO(child.communicate()[0].decode('utf-8')), 'hmmer3-text')
    except IOError:
        return print(f'\n[!] hmmer failed to run, is {in_file} DNA?')



def run_hmmbuild(in_file, id, out_file, threads):
    if which('hmmbuild') is None:
        exit('[!] hmmbuild not found')
    print("[>] Running hmmbuild...  ", end="", flush=True)
    cmd = ['hmmbuild', '--cpu', threads,
           '--informat', 'clustal', '--amino',
           '-n', id, out_file, '-']
    child = Popen(cmd, stdin=PIPE, stderr=DEVNULL, stdout=DEVNULL)
    child.stdin.write(in_file.encode())
    child.communicate()
    return out_file


def run_hmmpress(lib):
    if which('hmmpress') is None:
        exit('[!] hmmpress not found')
    print("[>] Running hmmpress...  ", end="", flush=True)
    Popen(["hmmpress", lib], stdout=DEVNULL)
    return print("done!")
