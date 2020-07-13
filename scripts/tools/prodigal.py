#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from subprocess import Popen, PIPE

def run_prodigal(in_file, lowqual):
    if lowqual is True:
        mode = 'meta'
    else:
        mode = 'single'
    print(f"[>] Running prodigal in {mode} mode...  ", end="", flush=True)
    cmd = ['prodigal', '-i', in_file, '-o', '/dev/null',
           '-a', '/dev/stdout', '-q', '-p', mode]
    child = Popen(cmd, stdout=PIPE)
    return child.communicate()[0].decode('utf-8')