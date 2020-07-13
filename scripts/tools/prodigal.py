#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from subprocess import Popen, PIPE


def run_prodigal(in_file, quality):
    print("Extracting proteins...")
    cmd = ['prodigal', '-i', in_file, '-o', '/dev/null',
           '-a', '/dev/stdout', '-q', '-p', quality]
    child = Popen(cmd, stdout=PIPE)
    return child.communicate()[0].decode('utf-8')