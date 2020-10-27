#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from gzip import open as gopen
from os import remove, rename

from sideroscanner.methods.fetch import fetch_url
from sideroscanner.tools.cdhit import run_cdhit
from sideroscanner.tools.blast import run_makeblastdb

def get_flankdb(flankpath, threads):
    Path(flankpath).mkdir(parents=True, exist_ok=True)
    print("[-] Preparing flanking virulence gene database")
    patric = fetch_url('ftp://ftp.patricbrc.org/specialty_genes/referenceDBs/PATRIC_VF.faa', None, flankpath + '/patric.faa')
    victors = fetch_url('http://www.phidias.us/victors/downloads/gen_downloads_protein.php', None, flankpath + '/victors.faa')
    vfdb = fetch_url('http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz', None, flankpath + '/vfdb.faa.gz')
    params = {'query': 'siderophore AND '
                       'taxonomy:"Bacteria [2]" AND '
                       'NOT receptor NOT partial NOT fragment', 'format': 'fasta'}
    bgcs = fetch_url('http://www.uniprot.org/uniprot/', params, flankpath + '/bgcs.faa')

    filenames = [patric, victors, vfdb, bgcs]
    with open(flankpath + '/flankdb', 'wb') as flank_file:
        for fname in filenames:
            if fname.endswith('.gz'):
                with gopen(fname) as infile:
                    for line in infile:
                        flank_file.write(line)
            else:
                with open(fname, 'rb') as infile:
                    for line in infile:
                        flank_file.write(line)
            remove(fname)
    run_makeblastdb(run_cdhit(flankpath + '/flankdb', 1, threads), 'prot', flankpath + '/flankdb')
    return rename(flankpath + '/flankdb_nr.faa', flankpath + '/flankdb')