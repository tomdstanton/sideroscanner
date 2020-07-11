#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__title__ = 'SideroScanner'
__version__ = '0.0.1'
__author__ = 'Tom Stanton'

import re
from sys import argv, stderr, exit, stdin
from os import getcwd, path, uname, cpu_count
from shutil import get_terminal_size, copyfileobj
from argparse import RawTextHelpFormatter, ArgumentParser
from datetime import datetime
from io import StringIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqIO import write, parse, to_dict
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import gzip
import tempfile

from scripts.config import hmmpath, plspath, mgepath, flankpath, furpath
from scripts.blast import run_blastn, run_blastp
from scripts.hmmer3 import run_hmmsearch, run_hmmpress, run_hmmscan
from scripts.prodigal import run_prodigal
from scripts.meme import run_mast
from scripts.fetch import fetch


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
    group.add_argument('--lib', metavar='hmm', type=str, default=hmmpath + '/iromps.hmm',
                       help='''| path/to/custom/HMM/(lib)rary.hmm
-----------------------------------------------''')
    group.add_argument('-v', action='store_true',
                       help='''| show version and exit
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''| show this help message and exit''')
    # if len(argv) == 1:
    #     parser.print_help(file=stderr)
    #     exit(1)
    return parser.parse_args()


def domain_filter(in_file, qual):
    q = run_hmmsearch(in_file, hmmpath+'/PF07715.hmm', threads)
    in_file_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    filter1 = ''
    for h in q.hit_keys:
        filter1 = filter1 + in_file_dict[h].format("fasta")
    print(f'Filtered {str(filter1.count(">"))} proteins with {q.description}')
    if qual == 'meta' or len(filter1) == 0:
        return filter1
    else:
        q = run_hmmsearch(filter1, hmmpath+'/PF00593.hmm', threads)
        filter2 = ''
        filter1_dict = to_dict(parse(StringIO(filter1), 'fasta'))
        for h in q.hit_keys:
            filter2 = filter2 + filter1_dict[h].format("fasta")
        print(f'Filtered {str(filter2.count(">"))} proteins with {q.description}')
        return filter2


def annotation(in_file, molecule):
    data = []
    for q in run_hmmscan(in_file, parse_args().lib, threads):
        if len(q.hits) > 0:
            data.append(f'{q.id.rsplit("_", 1)[0]},{q.id}'
                        f',{q.hits[0].id},{q.hits[0].description}'
                        f',{str(q.hits[0].bitscore)}')
        else:
            return None
    print(f'Annotated {len(data)} proteins')
    hmmscan_df = pd.DataFrame([sub.split(",") for sub in data],
                              columns=['contig', 'query', 'hit', 'description', 'score'])
    data = []
    for r in parse(StringIO(in_file), 'fasta'):
        data.append(f'{r.id},{str(len(r.seq))},'
                    f'{str(round((ProteinAnalysis(r._seq._data).molecular_weight()/1000),1))}')
    len_mass_df = pd.DataFrame([sub.split(",") for sub in data],
                               columns=['query', 'len', 'kDa'])
    if molecule == 'aa':
        return hmmscan_df.merge(len_mass_df, on='query', how='left')
    else:
        data = []
        for r in parse(StringIO(in_file), 'fasta'):
            data.append(','.join(r.description.replace(' ', '').split('#')[0:4]))
        bed_df = pd.DataFrame([sub.split(",") for sub in data],
                              columns=['query', 'start', 'end', 'str'])
        return hmmscan_df.merge(len_mass_df, on='query', how='left').merge(bed_df, on='query', how='left')


def range_subset(range1, range2):
    if not range1:
        return True
    if not range2:
        return False
    if len(range1) > 1 and range1.step % range2.step:
        return False
    return range1.start in range2 and range1[-1] in range2


def location(in_file, hits):
    print("Screening for plasmids...")
    plasmids = []
    percid = parse_args().l
    for q in run_blastn(in_file, plspath+'/plsdb.fna', 10, percid, threads):
        for h in q.hits:
            plasmids.append(f'{h.query_id}#{h.id}')
    print(f'{len(plasmids)} plasmid(s) found')
    print("Screening for MGEs...")
    mges = []
    for q in run_blastn(in_file, mgepath+'/mgedb', 0, percid, threads):
        for h in q.hits:
            mges.append(f'{h.query_id}#{h.id.split("|")[2]}#'
                        f'{str(h.hsps[0].query_range[0])}-'
                        f'{str(h.hsps[0].query_range[1])}')
    print("%i MGE(s) found"%len(mges))
    for row in hits.itertuples():
        for line in mges:
            start = int(line.split('#')[2].split('-')[0])
            end = int(line.split('#')[2].split('-')[1])
            if range_subset(range(int(row.start), int(row.end)), range(start, end)):
                hits.at[row.Index, 'plasmid/mge'] = line.split('#')[1]
    for row in hits.itertuples():
        for line in plasmids:
            if row.contig == line.split('#')[0]:
                hits.at[row.Index, 'plasmid/mge'] = line.split('#')[1]
    return hits


def flank_screen(in_file, hits):
    cds = parse_args().f
    cds_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    queries = ''
    for h in hits['query'].tolist():
        pos = int(h.split('_')[-1])
        acc = f'{h.rsplit("_", 1)[0]}_'
        for i in range(cds):
            i = i + 1
            try:
                u = cds_dict[acc + str(pos - i)]
                u.id = h
                u.description = f'up_{str(i)}'
                queries = queries + u.format("fasta")
            except:
                continue
            try:
                d = cds_dict[acc + str(pos + i)]
                d.id = h
                d.description = f'down_{str(i)}'
                queries = queries + d.format("fasta")
            except:
                continue
    up = []
    down = []
    print(f'Screening {cds} up/downstream CDS')
    for q in run_blastp(queries, flankpath+'/flankdb', '1e-130', threads):
        if len(q.hits) > 0:
            gene_name = re.sub("\s*\[[^[]*\]$", '', q.hsps[0].hit_description).strip()
            if gene_name.startswith('('):
                gene_name = gene_name.split('(', 1)[1].split(')')[0]
            if 'up' in q.hsps[0].query_description:
                up.append(f'{q.hsps[0].query_id}#{gene_name}')
            if 'down' in q.hsps[0].query_description:
                down.append(f'{q.hsps[0].query_id}#{gene_name}')
    u_df = pd.DataFrame([sub.split("#") for sub in up], columns=['query', 'upstream'])
    d_df = pd.DataFrame([sub.split("#") for sub in down], columns=['query', 'downstream'])
    u_df = u_df.groupby(['query'])['upstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    d_df = d_df.groupby(['query'])['downstream'].apply(lambda x: ' '.join(x.astype(str))).reset_index()
    return hits.merge(u_df, on='query', how='left').merge(d_df, on='query', how='left')

def tfbs_screen(in_file, hits):
    length = parse_args().b
    prom_dict = to_dict(parse(in_file, 'fasta'))
    queries = ''
    for row in hits.itertuples():
        if row.str == '1':
            f_start = int(row.start) - length
            f_end = int(row.start)
        elif row.str == '-1':
            f_start = int(row.end)
            f_end = int(row.end) + length
        try:
            x = prom_dict[row.contig][f_start:f_end]
            x.id = f'{row.query}:{str(f_start)}'
            queries = queries + x.format("fasta")
        except:
            continue
    print(f'Screening {length}bp upstream of hits for Fur binding sites')
    mast_out = run_mast(queries, furpath+'/fur.meme').rstrip()
    results = ''
    for line in (line for line in mast_out.split('\n') if not line.startswith('#')):
        results = results + line
    if len(results) == 0:
        print("No binding sites found")
        return hits
    else:
        mast = []
        for line in results.split('\n'):
            query = line.split(' ')[0].split(':')[0]
            fur_start = int(line.split(' ')[0].split(':')[1]) + int(line.split(' ')[4])
            fur_end = int(line.split(' ')[0].split(':')[1]) + int(line.split(' ')[5])
            pval = str(line.split(' ')[8])
            bs = prom_dict[query.rsplit('_',1)[0]][fur_start:fur_end].seq._data
            if line.split(' ')[1] == '-1':
                bs = str(Seq(bs).reverse_complement())
            mast.append(f'{query}#{str(fur_start)}#{str(fur_end)}#{pval}#{bs}')
        print(f'Putative Fur binding sites found for {len(mast)} hit(s)')
        df = pd.DataFrame([sub.split("#") for sub in mast],
                          columns=['query', 'fur_start',
                                   'fur_end', 'p_value',
                                   'fur_box'])
        return hits.merge(df, on='query', how='left')


def export_proteins(in_file, hits, seq):
    to_write = []
    in_file_dict = to_dict(parse(StringIO(in_file), 'fasta'))
    for row in hits.itertuples():
        r = in_file_dict[row.query]
        r.id = row.query
        r.description = f'{row.description} {row.hit}'
        to_write.append(r)
    if parse_args().e == 'seq':
        filename = f'{path.splitext(seq)[0]}_sideroscanner.faa'
    else:
        filename = parse_args().e
    write(to_write, filename, "fasta")
    print(f'Proteins written to: {filename}')


# def export_gff():
#     out_file = "your_file.gff"
#
#     qualifiers = {"source": "prediction", "score": 10.0, "other": ["Some", "annotations"],
#                   "ID": "gene1"}
#     sub_qualifiers = {"source": "prediction"}
#     top_feature = SeqFeature(FeatureLocation(0, 20), type="gene", strand=1,
#                              qualifiers=qualifiers)
#
#     rec.features = [top_feature]


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
            print(f"{seq} is not a file, skipping...")
            continue

        elif path.getsize(seq) <= 10:
            print(f"{seq} is too small, skipping...")
            continue

        if seq.endswith(".gz"):
            tmp = tempfile.NamedTemporaryFile()
            copyfileobj(gzip.open(seq), tmp)
            seq = tmp.name

        # Check molecule type and restrictions
        try:
            with open(seq, "r") as test:
                header = test.readline()
                firstseq = test.readline()
            if header.startswith(">"):
                pass
            else:
                print(f"{seq} is not a fasta file, skipping...")
                continue
        except:
            print(f"{seq} invalid file, skipping...")
            continue

        if "E" in firstseq:
            molecule = 'aa'
            print(f'Protein input detected: {name}')
        else:
            molecule = 'dna'
            print(f'DNA input detected: {name}')

        # If DNA input, get proteins
        quality = 'single'
        if molecule == 'dna':
            small_contig_list = []
            for r in parse(seq, "fasta"):
                if len(r.seq) < 10000:
                    small_contig_list.append(r.id)

            if len(small_contig_list) > 50:
                print(f'{str(len(small_contig_list))} contigs are < 10kbp')
                quality = 'meta'

            if len([1 for line in open(seq) if line.startswith(">")]) >= 1000:
                print(f'{name} has >= 1000 contigs')
                quality = 'meta'

            if parse_args().lowqual is not None:
                if len(parse_args().lowqual) == 0:
                    quality = 'meta'
                if seq in parse_args().lowqual:
                    quality = 'meta'

            if quality == 'meta':
                print('Switching Prodigal to meta mode')

            proteins = run_prodigal(seq, quality).replace('*', '')
            print(f"{proteins.count('>')} proteins extracted")

            if proteins.count('>') > 7500:
                print('Too many CDS, switching to lowqual mode')
                print('Re-running Prodigal in meta mode')
                quality = 'meta'
                proteins = run_prodigal(seq, quality).replace('*', '')
                print(f"{proteins.count('>')} proteins extracted") 

        else:
            proteins = ''
            for r in parse(seq, 'fasta'):
                proteins = proteins + r.format("fasta").replace('*', '')
            print(f"{proteins.count('>')} total protein queries")

        if proteins.count('>') == 0:
            continue

        # Filter
        iromps = domain_filter(proteins, quality)

        if len(iromps) == 0:
            print('No significant hits')
            continue
        # Annotate
        hits = annotation(iromps, molecule)
        if hits is None:
            print('No significant hits')
            continue

        if parse_args().f is not None:
            if molecule == 'aa':
                print('Cannot screen flanking CDS in protein fasta')
            else:
                hits = flank_screen(proteins, hits)

        if parse_args().l is not None:
            if molecule == 'aa':
                print('Cannot screen plasmids in protein fasta')
            else:
                hits = location(seq, hits)

        if parse_args().b is not None:
            if molecule == 'aa':
                print('Cannot screen Fur binding sites in protein fasta')
            else:
                hits = tfbs_screen(seq, hits)

        if parse_args().e is not None:
            export_proteins(iromps, hits, seq)

        # Tidying up
        hits = hits.drop(columns=['contig']).fillna("-")

        if parse_args().o is None:
            pd.set_option('display.max_colwidth', 10)
            print(hits.to_markdown(showindex=False))
        else:
            hits.insert(0, 'sample', name)
            out_df = out_df.append(hits)

    if parse_args().o is not None:
        out_df.to_csv(parse_args().o, index=False)
        print(f"Written results to: {parse_args().o}")


if __name__ == "__main__":
    if parse_args().t > cpu_count():
        print(f'Number of threads exceeds available CPUs, will use: {cpu_count()}')
        threads = str(cpu_count())
    else:
        threads = str(parse_args().t)
    main()
