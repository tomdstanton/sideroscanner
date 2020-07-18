#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez

def efetch(acc):
    print(f"[>] Fetching {acc} from NCBI...  ", end="", flush=True)
    Entrez.email = "tomdstanton@gmail.com"
    print('done!')
    return Entrez.efetch(
        db='protein', id=acc, rettype="fasta",
        retmode="text").read().rstrip('\n')