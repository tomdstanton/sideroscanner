#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from shutil import which

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None

# Check if tools are installed
    for d in ['hmmscan', 'hmmsearch', 'blastp', 'blastn', 'prodigal']:
        if not is_tool(d):
            print('ERROR: ' + d + ' not found')
            sys.exit(1)