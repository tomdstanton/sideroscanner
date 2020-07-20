#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from subprocess import Popen, PIPE
from shutil import which
from sys import exit

def run_prodigal(in_file, lowqual):
    if which('prodigal') is None:
        exit('[!] prodigal not found')

    if lowqual is True:
        mode = 'meta'
    else:
        mode = 'single'
    print(f"[>] Running prodigal in {mode} mode...  ", end="", flush=True)
    cmd = ['prodigal', '-i', in_file, '-o', '/dev/null',
           '-a', '/dev/stdout', '-q', '-p', mode]
    child = Popen(cmd, stdout=PIPE)
    return child.communicate()[0].decode('utf-8')