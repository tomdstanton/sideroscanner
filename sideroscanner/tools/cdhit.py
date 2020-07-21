#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from subprocess import run, DEVNULL
from os import path, remove
from sys import exit
from shutil import which

def run_cdhit(in_file, perc):
    if which('cd-hit') is None:
        exit('[!] cd-hit not found')
    out_file = path.splitext(in_file)[0]+'_nr.faa'
    print("[>] Running cd-hit...    ", end="", flush=True)
    cmd = ['cd-hit', '-i', in_file, '-o', out_file, '-c', str(perc), '-T', '0']
    try:
        run(cmd, stdout=DEVNULL)
        print("done!")
        remove(out_file + '.clstr')
        remove(in_file)
        return out_file
    except:
        raise SystemExit()