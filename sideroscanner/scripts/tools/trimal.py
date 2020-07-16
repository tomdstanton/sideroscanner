#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from subprocess import Popen, PIPE
from shutil import which
from sys import exit

def run_trimal(in_file):
    if which('trimal') is None:
        exit('[!] trimal not found')
    print('[>] Running trimal...    ', end="", flush=True)
    cmd = ['trimal', '-in', in_file, '-automated1', '-clustal']
    child = Popen(cmd, stdout=PIPE)
    print('done!')
    return child.communicate()[0].decode('utf-8')