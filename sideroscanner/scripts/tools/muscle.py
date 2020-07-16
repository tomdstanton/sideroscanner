#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from subprocess import Popen, PIPE, DEVNULL
from shutil import which
from sys import exit

def run_muscle(in_file, out_file):
    if which('muscle') is None:
        exit('[!] muscle not found')
    else:
        print("[>] Running muscle...  ", end="", flush=True)
        cmd = ['muscle', '-out', out_file]
        child = Popen(cmd, stderr=DEVNULL, stdin=PIPE)
        child.stdin.write(in_file.encode())
        child.communicate()
        print('done!')
        return out_file