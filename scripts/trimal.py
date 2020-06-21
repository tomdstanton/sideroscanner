#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from subprocess import Popen, PIPE

def run_trimal(in_file):
    cmd = ['trimal', '-in', in_file, '-automated1', '-clustal']
    child = Popen(cmd, stdout=PIPE)
    return child.communicate()[0].decode('utf-8')