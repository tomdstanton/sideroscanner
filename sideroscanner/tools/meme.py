#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from subprocess import Popen, PIPE, DEVNULL
from shutil import which
from sys import exit

def run_mast(in_file, pwm):
    if which('mast') is None:
        exit('[!] mast not found')
    print("[>] Running mast...  ", end="", flush=True)
    cmd = ['mast', '-hit_list', '-best', pwm, '-']
    child = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=DEVNULL)
    child.stdin.write(in_file.encode())
    return child.communicate()[0].decode('utf-8')