#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from subprocess import Popen, PIPE
from shutil import which
from sys import exit
from tempfile import SpooledTemporaryFile as tempfile

def run_prodigal(in_file, lowqual):
    if which('prodigal') is None:
        exit('[!] prodigal not found')

    if lowqual is True:
        mode = 'meta'
    else:
        mode = 'single'
    print(f"[>] Running prodigal in {mode} mode...  ", end="", flush=True)
    cmd = ['prodigal', '-i', '/dev/stdin', '-o', '/dev/null',
           '-a', '/dev/stdout', '-q', '-p', mode]
    f = tempfile()
    f.write(in_file.encode())
    f.seek(0)
    child = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=f)
    return child.communicate()[0].decode('utf-8')
#f.close()