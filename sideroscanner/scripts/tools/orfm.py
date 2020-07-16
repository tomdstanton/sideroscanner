#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from subprocess import Popen, PIPE
from shutil import which
from sys import exit

def run_orfm(in_file):
    if which('orfm') is None:
        exit('[!] orfm not found')
    print(f"[>] Running orfm...  ", end="", flush=True)
    cmd = ['orfm', in_file, '-c', '11']
    child = Popen(cmd, stdout=PIPE)
    return child.communicate()[0].decode('utf-8')