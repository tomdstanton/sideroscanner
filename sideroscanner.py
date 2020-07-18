#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import tempfile
from argparse import RawTextHelpFormatter, ArgumentParser
from datetime import datetime
from os import path, uname, cpu_count, stat
from os.path import dirname, abspath
from shutil import get_terminal_size, copyfileobj
from sys import exit, stdin, stderr, argv

import pandas as pd
from Bio.SeqIO import parse

from scripts.annotation import domain_filter, annotate
from scripts.export import export_proteins
from scripts.fetch import fetch_url
from scripts.flanks import flank_screen
from scripts.location import plasmid_mge_screen
from scripts.tfbs import tfbs_screen
from scripts.tools.hmmer3 import run_hmmpress
from scripts.tools.orfm import run_orfm
from scripts.tools.prodigal import run_prodigal

dbpath = abspath(dirname(argv[0]))+'/data'
mgepath = dbpath+'/mgedb'
plspath = dbpath+'/plsdb'
flankpath = dbpath+'/flankdb'
furpath =  dbpath+'/furdb'
iromppath = dbpath+'/irompdb'
blastdb = iromppath+'/irompdb'
lib = iromppath+"/iromps.hmm"

def parse_args():
    parser = ArgumentParser(add_help=False,
                                     formatter_class=RawTextHelpFormatter,
                                     usage="sideroscanner.py -i query.fasta",
                                     description='''
SideroScannner: a tool for annotating IROMPs in bacteria
-----------------------------------------------------------------
Please cite: Stanton et al, 2020''')
    group = parser.add_argument_group("Options")

    group.add_argument('-i', metavar='-', nargs='*', type=str,
                       help='''path/to/input ['-' for STDIN]
-----------------------------------------------''')
    group.add_argument('-o', metavar='-', nargs='?', type=str,
                       const=f'sideroscanner_{datetime.today().strftime("%d%m%y_%H%M")}.csv',
                       help='''output results.csv instead of STDOUT
[optional: path/to/output/file]
[default: sideroscanner_DDMMYY_hhmm.csv]
-----------------------------------------------''')
    group.add_argument('-l', metavar='int', nargs='?', type=int, const=90,
                       help='''genomic location screen
[optional: blastn percid] [default: 90]
-----------------------------------------------''')
    group.add_argument('-f', metavar='int', nargs='?', type=int, const=3,
                       help='''flanking CDS screen
[optional: number of up/downstream CDS]
[default: 3]
-----------------------------------------------''')
    group.add_argument('-b', metavar='int', nargs='?', type=int, const=200,
                       help='''Fur binding site screen
[optional: length of promoter region]
[default: 200]
-----------------------------------------------''')
    group.add_argument('-e', metavar='-', nargs='?', type=str,
                       const='',
                       help='''export annotated proteins
[optional: path/to/export/fasta]
[default: filename_sideroscanner.faa]
-----------------------------------------------''')
    group.add_argument('-t', metavar='int', type=int, default=cpu_count(),
                       help=f'''number of threads [default: {cpu_count()}]
-----------------------------------------------''')
    group.add_argument('--lowqual', metavar='-', nargs='?', type=str, const='',
                       help=''''meta' CDS prediction AND single domain filter
[optional: path/to/draft/genome]
[default: all inputs]
-----------------------------------------------''')
    group.add_argument('--nofilter', metavar='-', nargs='?', type=str, const='',
                   help='''turn off domain filter
[optional: path/to/input/fasta]
[default: all inputs]
-----------------------------------------------''')
    group.add_argument('--lib', metavar='hmm', type=str, default=lib,
                       help='''path/to/custom/HMM/library.hmm
-----------------------------------------------''')
    group.add_argument('-v', action='store_true',
                       help='''show version and exit
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''show this help message and exit''')
    if len(argv) == 1:
        parser.print_help(file=stderr)
        exit(1)
    return parser.parse_args()


def main():
    if parse_args().v is True:
        exit(f'{__version__}')

    if len(parse_args().i) == 0:
        exit('[!] Please enter a valid input with [-i]')

    # Some quick setup checks
    # Check database folder
    if not path.isdir(iromppath):
        exit(f'[!] {iromppath} is not an existing directory')
    else:
        plug = iromppath + '/PF07715.hmm'
        if not path.isfile(plug):
            print(f'[!] {plug} not found, will download')
            fetch_url('https://pfam.xfam.org/family/PF07715/hmm', None, plug)
        receptor = iromppath + '/PF00593.hmm'
        if not path.isfile(receptor):
            print(f'[!] {receptor} not found, will download')
            fetch_url('https://pfam.xfam.org/family/PF00593/hmm', None, receptor)

    # Check library HMM
    lib = parse_args().lib
    if not path.isfile(lib):
        exit(f'[!] {lib} is not a valid file \n'
             f'[!] Try running: ./sideroscanner-buildhmms.py')
    else:
        for i in [".h3i", ".h3p", ".h3f", ".h3m"]:
            if not path.isfile(f'{lib}{i}'):
                print(f'[!] {lib} is not pressed')
                run_hmmpress(lib)
                break

    if parse_args().f is not None:
        if not path.isfile(flankpath+'/flankdb.psq'):
            exit(f'[!] {flankpath}/flankdb.psq not found\n'
                 f'[!] Try running: ./sideroscanner-builddbs.py -db flankdb')

    if parse_args().l is not None:
        if not path.isfile(plspath+'/plsdb.fna.nsq') and not path.isfile(mgepath+'/mgedb.nsq'):
            exit(f'[!] {plspath}/plsdb.fna.nsq and {mgepath}/mgedb.nsq not found\n'
                 f'[!] Try running: ./sideroscanner-builddbs.py -db plsdb mgedb')

    if parse_args().b is not None:
        if not path.isfile(furpath+'/fur.meme'):
            exit(f'[!] {furpath}/fur.meme not found')

    # Sanity checks passed, assuming good to go
    print('-' * int(get_terminal_size()[0]))
    print(f'{__name__} {__version__}')
    print(f'Your system is {uname()[0]}')
    if 'Linux' not in uname()[0]:
        print(f'[!] Warning: scripts has not been tested on {uname()[0]}')
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
        if parse_args().lowqual is not None:
            if len(parse_args().lowqual) == 0 or seq in parse_args().lowqual:
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
        if parse_args().nofilter is not None:
            if len(parse_args().nofilter) == 0 or seq in parse_args().nofilter:
                iromps = proteins
            else: iromps = domain_filter(proteins, iromppath, lowqual, threads)
        else: iromps = domain_filter(proteins, iromppath, lowqual, threads)

        if len(iromps) == 0:
            print('[!] No significant hits')
            continue
        # Annotate
        hits = annotate(iromps, input_type, parse_args().lib, threads)
        if hits is None:
            print('[!] No significant hits')
            continue

        if parse_args().f is not None:
            if input_type is not 'genome':
                print(f'[!] Cannot screen flanking CDS in {input_type} input file')
            else:
                hits = flank_screen(proteins, hits, flankpath, parse_args().f, threads)

        if parse_args().l is not None:
            if input_type is not 'genome':
                print(f'[!] Cannot screen plasmids/MGEs in {input_type} input file')
            else:
                hits = plasmid_mge_screen(seq, plspath, mgepath, hits, parse_args().l, threads)

        if parse_args().b is not None:
            if input_type is not 'genome':
                print(f'[!] Cannot screen Fur binding sites in {input_type} input file')
            else:
                hits = tfbs_screen(seq, furpath, hits, parse_args().b)

        if parse_args().e is not None:
            export_proteins(iromps, parse_args().e, name, hits)

        # Tidying up
        hits = hits.drop(columns=['contig']).fillna("-")

        if parse_args().o is None:
            print(hits.to_markdown(showindex=False))

        else:
            hits.insert(0, 'sample', name)
            out_df = out_df.append(hits)

    if parse_args().o is not None and not out_df.empty:
        out_df.to_csv(parse_args().o, index=False)
        print(f"[-] Written results to: {parse_args().o}")

if __name__ == "__main__":
    if parse_args().t > cpu_count():
        print(f'[!] Number of threads exceeds available CPUs, will use: {cpu_count()}')
        threads = str(cpu_count())
    else:
        threads = str(parse_args().t)
    main()
