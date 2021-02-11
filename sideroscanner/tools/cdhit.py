#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from subprocess import Popen, PIPE, DEVNULL
from sys import exit
from shutil import which
from tempfile import SpooledTemporaryFile as tempfile

def run_cdhit(in_file, perc):
    if which('cd-hit') is None:
        exit('[!] cd-hit not found')

    print(f"[>] Clustering at {perc}%...    ", end="", flush=True)
    cmd = ['cd-hit', '-i', '/dev/stdin', '-o', '/dev/stdout', '-c', str(perc/100), '-M', '0']

    try:
        f = tempfile()
        f.write(in_file.encode())
        f.seek(0)
        child = Popen(cmd, stderr=DEVNULL, stdin=f, stdout=PIPE)
        proteins = child.communicate()[0].decode('utf-8')
        num_prot = proteins.count('>')
        print(f"{num_prot} proteins remain")
        return f">{proteins.split('>',1)[1].split('     ', -1)[0].strip()}"
    except:
        exit('CD-HIT ERROR')