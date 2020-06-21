#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from subprocess import Popen, PIPE

def run_muscle(in_file, out_file):
    cmd = ['muscle', '-out', out_file]
    child = Popen(cmd, stdout=PIPE, stdin=PIPE)
    child.stdin.write(in_file.encode())
    child.communicate()
    return out_file