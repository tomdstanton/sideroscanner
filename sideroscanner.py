#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import gzip
import tempfile
from argparse import RawTextHelpFormatter, ArgumentParser
from datetime import datetime
from os import getcwd, path, uname, cpu_count, stat
from shutil import get_terminal_size, copyfileobj
from sys import exit, stdin, stderr, argv

import pandas as pd
from Bio.Seq import translate
from Bio.SeqIO import parse

import warnings
from Bio import BiopythonWarning

from scripts import tfbs_screen, location, flank_screen, domain_filter, annotation
from scripts.export import export_proteins
from scripts.config import hmmpath
from scripts.fetch import fetch
from scripts.tools.hmmer3 import run_hmmpress
from scripts.tools.prodigal import run_prodigal
from scripts.tools.orfm import run_orfm


def parse_args():
    parser = ArgumentParser(add_help=False,
                                     formatter_class=RawTextHelpFormatter,
                                     usage="sideroscanner.py -i query.fasta",
                                     description='''
        SideroScannner: a tool for annotating IROMPs in bacteria
        ========================================================
        Please cite: Stanton et al, 2020''')
    group = parser.add_argument_group("Options")

    group.add_argument('-i', metavar='-', nargs='*', type=str,
                       help='''| path/to/(i)nput/fasta ['-' for STDIN]
-----------------------------------------------''')
    group.add_argument('-o', metavar='-', nargs='?', type=str,
                       const=f'sideroscanner_{datetime.today().strftime("%d%m%y_%H%M")}.csv',
                       help=f'''| (o)utput file.csv instead of STDOUT
    [optional: path/to/(o)utput/file]
    [default: {getcwd()}/sideroscanner_DDMMYY_hhmm.csv]
-----------------------------------------------''')
    group.add_argument('-l', metavar='int', nargs='?', type=int, const=90,
                       help='''| determine genomic (l)ocation of hits
    [optional: blastn percid]
    [default: 90]
-----------------------------------------------''')
    group.add_argument('-f', metavar='int', nargs='?', type=int, const=3,
                       help='''| (f)lanking CDS screen
    [optional: number of up/downstream CDS]
    [default: 3]
-----------------------------------------------''')
    group.add_argument('-b', metavar='int', nargs='?', type=int, const=200,
                       help='''| Fur (b)inding site screen
    [optional: length of promoter region]
    [default: 200]
-----------------------------------------------''')
    group.add_argument('-e', metavar='-', nargs='?', type=str,
                       const='seq',
                       help='''| (e)xport annotated proteins
    [optional: path/to/export/fasta]
    [default: path/to/'input'_sideroscanner.faa]
-----------------------------------------------''')
    group.add_argument('-t', metavar='int', type=int, default=cpu_count(),
                       help=f'''| number of (t)hreads to use
    [default: {cpu_count()}]
-----------------------------------------------''')
    group.add_argument('--lowqual', metavar='-', nargs='?', type=str, const='',
                       help='''| (1) 'meta' CDS prediction AND (2) filters with plug domain only
    [optional: path/to/(low)/(qual)ity/input/fasta]
    [default: all inputs]
-----------------------------------------------''')
    group.add_argument('--nofilter', metavar='-', nargs='?', type=str, const='',
                   help='''| Turn off domain filter
    [optional: path/to/input/fasta]
    [default: all inputs]
-----------------------------------------------''')
    group.add_argument('--lib', metavar='hmm', type=str, default=hmmpath + '/iromps.hmm',
                       help='''| path/to/custom/HMM/(lib)rary.hmm
-----------------------------------------------''')
    group.add_argument('-v', action='store_true',
                       help='''| show version and exit
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''| show this help message and exit''')
    if len(argv) == 1:
        parser.print_help(file=stderr)
        exit(1)
    return parser.parse_args()


def main():
    if parse_args().v is True:
        print(f'{__title__} {__version__}')

    # Some quick setup checks
    # Check database folder
    if not path.isdir(hmmpath):
        print(f'{hmmpath} is not an existing directory')
        exit(1)
    else:
        plug = hmmpath+'/PF07715.hmm'
        if not path.isfile(plug):
            print(f'{plug} not found, will download')
            fetch('https://pfam.xfam.org/family/PF07715/hmm', None, plug)
        receptor = hmmpath+'/PF00593.hmm'
        if not path.isfile(receptor):
            print(f'{receptor} not found, will download')
            fetch('https://pfam.xfam.org/family/PF00593/hmm', None, receptor)

    # Check library HMM
    lib = parse_args().lib
    if not path.isfile(lib):
        print(f'{lib} is not a valid file')
        print('You can download the current IROMP HMM library from:')
        print('https://github.com/tomdstanton/sideroscanner/blob/master/databases/HMMs/iromps.hmm')
        exit(1)
    else:
        for i in [".h3i", ".h3p", ".h3f", ".h3m"]:
            if not path.isfile(f'{lib}{i}'):
                print(f'{lib} is not pressed, allow me...')
                run_hmmpress(lib)
                break

    # Sanity checks passed, assuming good to go
    print('-' * int(get_terminal_size()[0]))
    print(f'{__title__} {__version__}')
    print(f'Your system is {uname()[0]}')
    if 'Linux' not in uname()[0]:
        print(f'Warning: sideroscanner has not been tested on {uname()[0]}')
    print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))
    print(f'Using {threads} threads...')

    if parse_args().o is not None:
        out_df = pd.DataFrame()

    # Loop over each input file
    for seq in parse_args().i:
        print('-' * int(get_terminal_size()[0]))
        if seq == '-':
            name = 'stdin'
            tmp = tempfile.NamedTemporaryFile()
            with open(tmp.name, 'w') as f:
                for line in stdin:
                    f.write(line)
            seq = tmp.name
        else:
            name = path.splitext(path.basename(seq))[0]

        # Check valid inputs
        if not path.isfile(seq):
            print(f"[!] {seq} is not a file, skipping...")
            continue

        elif stat(seq).st_size <= 10:
            print(f"[!] {seq} is too small, skipping...")
            continue

        if seq.endswith(".gz"):
            tmp = tempfile.NamedTemporaryFile()
            copyfileobj(gzip.open(seq), tmp)
            original_seq = seq
            seq = tmp.name
            name = path.splitext(name)[0]

        # Check input_type type and restrictions
        try:
            with open(seq, "r") as test:
                header = test.readline()
                firstseq = test.readline()
            if header.startswith(">") or header.startswith("@"):
                pass
            else:
                print(f"[!] {seq} is not a fasta file, skipping...")
                continue
        except:
            print(f"[!] {seq} invalid file, skipping...")
            continue

        lowqual = False
        if parse_args().lowqual is not None and seq in parse_args().lowqual:
            lowqual = True

        if "E" in firstseq:
            input_type = 'protein'
            print(f'[-] Guessing {input_type} input: \033[4m{name}\033[0m')
            proteins = ''
            for r in parse(seq, 'fasta'):
                proteins = proteins + r.format("fasta").replace('*', '')
        elif '@' in header:
            input_type = 'fastq'
            print(f'[-] Guessing {input_type} input: \033[4m{name}\033[0m')
            proteins = run_orfm(seq)

        else:
            lengths = []
            for r in parse(seq, 'fasta'):
                lengths.append(len(r.seq) < 10000)

            if not all(lengths) is True and stat(seq).st_size >= 531000:
                input_type = 'genome'
                print(f'[-] Guessing {input_type} input: \033[4m{name}\033[0m')
                if lengths.count(True) > 50 or len(lengths) > 1000:
                    lowqual = True

                proteins = run_prodigal(seq, lowqual).replace('*', '')
                if proteins.count('>') > 7500:
                    print('\n[!] Too many CDS, switching to low quality mode')
                    lowqual = True
                    proteins = run_prodigal(seq, lowqual).replace('*', '')
            else:
                input_type = 'gene'
                print(f'[-] Guessing {input_type} input: \033[4m{name}\033[0m')
                proteins = run_orfm(seq)

            print(f"{proteins.count('>')} total protein queries")

        if proteins.count('>') == 0:
            continue

        # Filter
        if parse_args().nofilter is not None and seq in parse_args().nofilter:
            iromps = proteins
        else:
            iromps = domain_filter.domain_filter(proteins, lowqual, threads)

        if len(iromps) == 0:
            print('[!] No significant hits')
            continue
        # Annotate
        hits = annotation.annotation(iromps, input_type, parse_args().lib, threads)
        if hits is None:
            print('[!] No significant hits')
            continue

        if parse_args().f is not None:
            if input_type is not 'genome':
                print(f'[!] Cannot screen flanking CDS in {input_type} input file')
            else:
                hits = flank_screen.flank_screen(proteins, hits, parse_args().f, threads)

        if parse_args().l is not None:
            if input_type is not 'genome':
                print(f'[!] Cannot screen plasmids/MGEs in {input_type} input file')
            else:
                hits = location.location(seq, hits, parse_args().l, threads)

        if parse_args().b is not None:
            if input_type is not 'genome':
                print(f'[!] Cannot screen Fur binding sites in {input_type} input file')
            else:
                hits = tfbs_screen.tfbs_screen(seq, hits, parse_args().b)

        if parse_args().e is not None:
            if 'original_seq' in locals():
                seq = path.splitext(original_seq)[0]
            export_proteins(iromps, hits, seq)

        # Tidying up
        hits = hits.drop(columns=['contig']).fillna("-")

        if parse_args().o is None:
            pd.set_option('display.max_colwidth', 10)
            print(hits.to_markdown(showindex=False))
        else:
            hits.insert(0, 'sample', name)
            out_df = out_df.append(hits)

    if parse_args().o is not None and not out_df.empty:
        out_df.to_csv(parse_args().o, index=False)
        print(f"[-] Written results to: {parse_args().o}")


if __name__ == "__main__":
    if parse_args().t > cpu_count():
        print(f'Number of threads exceeds available CPUs, will use: {cpu_count()}')
        threads = str(cpu_count())
    else:
        threads = str(parse_args().t)
    main()
