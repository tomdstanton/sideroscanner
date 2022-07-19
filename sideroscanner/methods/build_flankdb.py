#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from gzip import open as gopen
from os import remove
from Bio.SeqIO import parse, write, to_dict
from io import StringIO

from sideroscanner.methods.fetch import fetch_url
from sideroscanner.tools.blast import run_makeblastdb


def get_flankdb(flankpath):
    Path(flankpath).mkdir(parents=True, exist_ok=True)

    print("[-] Preparing flanking virulence gene database")

    patric = fetch_url(
        'ftp://ftp.patricbrc.org/specialty_genes/referenceDBs/PATRIC_VF.faa', None, flankpath + '/patric.faa')
    # victors = fetch_url(
    #     'http://www.phidias.us/victors/downloads/gen_downloads_protein.php', None, flankpath + '/victors.faa')
    vfdb = fetch_url(
        'http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz', None, flankpath + '/vfdb.faa.gz')

    params = {'query': 'siderophore AND '
                       'taxonomy_id:2 AND '
                       'NOT receptor NOT partial NOT fragment', 'format': 'fasta'}
    bgcs = fetch_url('http://www.uniprot.org/uniprot/', params, flankpath + '/bgcs.faa')

    filenames = [patric, victors, vfdb, bgcs]

    db = ''

    for fname in filenames:
        if fname.endswith('.gz'):
            with gopen(fname, 'rt') as infile:
                for line in infile:
                    db += line
        else:
            with open(fname, 'rt') as infile:
                for line in infile:
                    db += line
        remove(fname)

    d1 = db.count('>')
    print(f"[-] {d1} total proteins downloaded")
    accessions = set()
    db2 = ''
    for r in parse(StringIO(db), 'fasta'):
        if r.id not in accessions:
            accessions.add(r.id)
            db2 += r.format('fasta')
    d2 = db2.count('>')
    print(f"[-] Removed {d1 - d2} duplicate accessions")

    fasta_lines = db2.split('>')[0:]  # splits each sequence by header

    def remove_complete_duplicates(fasta_lines):
        print(f"[>] Removing redundancy...  ", end="", flush=True)
        outputlist, setofuniqsequence = [], set()
        for sequence in fasta_lines:
            if sequence not in setofuniqsequence:
                outputlist.append(sequence)
                setofuniqsequence.add(sequence)
        print(f"{len(outputlist)} proteins remaining")
        return outputlist

    with open(flankpath + '/flankdb', 'w')  as flank_file:
        flank_file.write('>'.join(remove_complete_duplicates(fasta_lines)))

    return run_makeblastdb(flankpath + '/flankdb', 'prot', flankpath + '/flankdb')
