#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .fetch import fetch_url
from shutil import disk_usage
from sys import exit
from os import remove

def get_fastadb(iromppath):
    if (disk_usage("/")[2] // (2**30)) < 1:
        exit("[!] Need at least 1Gb free space!")
    else:
        try:
            url = 'http://www.uniprot.org/uniprot/'
            params = {'query': 'taxonomy:"Bacteria [2]" AND '
                               #'keyword:"TonB box [KW-0798]" OR '
                               'family:"tonb-dependent receptor family" '
                               'NOT partial NOT fragment', 'format': 'fasta'}
            return fetch_url(url, params, iromppath + '/iromps.faa')
        
        except KeyboardInterrupt:
            remove(iromppath+'/iromps.faa')
            exit(f'[!] Keyboard interrupt, removing {iromppath}/iromps.faa')

