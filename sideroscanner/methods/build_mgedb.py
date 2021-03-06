#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from Bio.SeqIO import parse, write
from os import remove
from io import StringIO

from sideroscanner.methods.fetch import fetch_url
from sideroscanner.tools.blast import run_makeblastdb

def get_mgedb(mgepath):
    print("[-] Preparing mobile genetic element database")
    Path(mgepath).mkdir(parents=True, exist_ok=True)

    to_write = []

    for i in range(1, 15):
        file = fetch_url(f'https://raw.githubusercontent.com/katholt/'
                         f'Kleborate/master/kleborate/data/ICEKp_references'
                         f'/ICEKp{i}.embl',
                         None,
                         f'{mgepath}/icekp{i}')

        flist = open(file).readlines()
        parsing = False
        fasta = ''
        for line in flist:
            if line.startswith("//"):
                parsing = False
            if parsing:
                fasta += line.replace(" ", "").strip()
            if line.startswith("SQ"):
                parsing = True
        icekp = f'>ICEKp{i}\n' + ''.join([i for i in fasta if not i.isdigit()])
        for r in parse(StringIO(icekp), 'fasta'):
            to_write.append(r)
        remove(file)

    ice = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/ICE_seq_all.fas', None, mgepath + '/ice.fna')
    t4ss = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/T4SS-type_ICE_seq_all.fas', None, mgepath + '/t4ss.fna')
    aice = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/AICE_seq_all.fas', None, mgepath + '/aice.fna')
    ime = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/IME_seq_all.fas', None, mgepath + '/ime.fna')
    cime = fetch_url('https://db-mml.sjtu.edu.cn/ICEberg2/download/CIME_seq_all.fas', None, mgepath + '/cime.fna')

    filenames = [ice, t4ss, aice, ime, cime]
    accessions = ['ICEKp1']
    for f in filenames:
        for r in parse(f, 'fasta'):
            r.id = r.id.split('|')[2]
            r.id = r.id.replace('[', '_')
            r.id = r.id.replace(']', '')
            if r.id not in accessions:
                accessions.append(r.id)
                to_write.append(r)
        remove(f)




    write(to_write, mgepath + '/mgedb', "fasta")

    return run_makeblastdb(mgepath + '/mgedb','nucl', f'{mgepath}/mgedb')