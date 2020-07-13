#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from subprocess import Popen, PIPE

def run_orfm(in_file):
    print(f"[>] Running orfm...  ", end="", flush=True)
    cmd = ['orfm', in_file, '-c', '11']
    child = Popen(cmd, stdout=PIPE)
    return child.communicate()[0].decode('utf-8')