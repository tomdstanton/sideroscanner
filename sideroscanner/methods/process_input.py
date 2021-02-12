#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.SeqIO import parse
from io import StringIO
from sys import exit, stdin
from gzip import open as gzopen
from sideroscanner.tools.orfm import run_orfm
from sideroscanner.tools.prodigal import run_prodigal
from pathlib import Path

def check_seq(file):
    if file == '-':
        try:
            seq = stdin.read()
        except UnicodeDecodeError:
            exit("[!] Cannot read STDIN, if gzipped try: file | gunzip -c | sideroscanner")
    else:
        if not Path(file).exists():
            return print(f"[!] {file} does not exist, skipping...")
        try:
            with open(file, 'r') as fr:
                seq = fr.read()
        except IsADirectoryError:
            return print(f"[!] {file} is a directory, skipping...")

        except UnicodeDecodeError:
            try:
                with gzopen(file, "rt") as fr:
                    seq = fr.read()
            except OSError:
                return print(f"[!] Could not open {file}, skipping...")

    if len(seq) <= 10:
        return print(f"[!] {file} is too small, skipping...")

    if ">" in seq[0] or "@" in seq[0]:
        return seq
    else:
        return print(f"[!] {file} is not a fasta file, skipping...")


def proteins_from_seq(seq, qual_arg, mol_arg):
    proteins = ''
    if ">" in seq[0]:
        lowqual = False
        if mol_arg == 'gene':
            input_type = 'gene'
            print(f'[-] User specified {input_type} input')
            proteins = run_orfm(seq)

        # Check if protein
        if "E" in '\n'.join(seq.split('\n')[1:])[0:99] or mol_arg == 'protein':
            input_type = 'protein'
            if mol_arg == 'protein':
                print(f'[-] User specified {input_type} input, ', end="", flush=True)
            else:
                print(f'[?] Guessing {input_type} input, ', end="", flush=True)
            proteins = seq.replace('*', '')

        # Assuming DNA
        else:
            if qual_arg is not None:
                if len(qual_arg) == 0 or seq in qual_arg:
                    lowqual = True

            length = [len(r.seq) for r in parse(StringIO(seq), 'fasta')]
            genome_calc = (sum(length) / len(length)) > 1000 and len(seq) >= 531000

            if genome_calc is True or mol_arg == 'genome':
                input_type = 'genome'
                if mol_arg == 'genome':
                    print(f'[-] User specified {input_type} input')
                else:
                    print(f'[?] Guessing {input_type} input')

                if len([g < 10000 for g in length]) > 50 or len(length) > 1000:
                    lowqual = True

                proteins = run_prodigal(seq, lowqual).replace('*', '')

                if proteins.count('>') > 7500 and lowqual is False:
                    print('\n[!] Too many CDS, switching to low quality mode')
                    lowqual = True
                    proteins = run_prodigal(seq, lowqual).replace('*', '')

            else:
                input_type = 'gene'
                print(f'[?] Guessing {input_type} input')
                proteins = run_orfm(seq)

    elif '@' in '\n'.join(seq.split('\n')[1:])[0:99]:
        input_type = 'fastq'
        print(f'[?] Guessing {input_type} input'
              f'[!] finding ORFs in fastq files depends on read length')
        proteins = run_orfm(seq)

    num_prot = proteins.count('>')
    if num_prot > 0:
        print(f"{num_prot} protein queries")
        return proteins, input_type, lowqual
    else:
        return print(f"\n[!] No proteins to scan, skipping...")


