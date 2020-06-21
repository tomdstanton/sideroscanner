#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

from subprocess import run
import os

def run_cdhit(in_file, perc):
    out_file = os.path.splitext(in_file)[0]+'_nr.faa'
    print('Clutsering ' + os.path.basename(in_file) + ' at '+str(perc*100)+'%...')
    try:
        run(['cd-hit', '-i', in_file, '-o', out_file, '-c', str(perc), '-T', '0'])
        os.remove(out_file + '.clstr')
        os.remove(in_file)
    except:
        raise SystemExit()
    return out_file