#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from subprocess import Popen, PIPE, DEVNULL

def run_mast(in_file, pwm):
    print("Running MAST...")
    cmd = ['mast', '-hit_list', '-best', pwm, '-']
    child = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=DEVNULL)
    child.stdin.write(in_file.encode())
    return child.communicate()[0].decode('utf-8')