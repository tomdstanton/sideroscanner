import sys
import os
import argparse
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd
import re
from argparse import RawTextHelpFormatter
from Bio import SeqIO, Entrez, SearchIO
from build_dbs import fetch, cluster, make_blast_db



pathname = os.path.dirname(sys.argv[0])
full_path = os.path.abspath(pathname)

def parse_args():
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=RawTextHelpFormatter,
                                     usage="./iromp_hmm_builder.py",
                                     description='''
        sideroscannner: a tool for annotating IROMPs in bacteria
        ========================================================''')
    group = parser.add_argument_group("Options")

    group.add_argument('-i', metavar='csv', type=str,
                       default=full_path + '/databases/iromps.csv',
                        help='''| path/to/iromp/csv
    [default: ''' + full_path + '''/databases/iromps.csv]
-----------------------------------------------''')
    group.add_argument('--dbpath', metavar='path', type=str, default=full_path + '/databases/',
                        help='''| path/to/db/
    [default: ''' + full_path + '''/databases/]
-----------------------------------------------''')
    group.add_argument("-h", action="help", help='''| show this help message and exit''')
    return parser.parse_args()

dbpath = parse_args().dbpath


def get_fastadb():
    url = 'http://www.uniprot.org/uniprot/'
    params = {'query':'taxonomy:"Bacteria [2]" AND ' \
        'family:"tonb-dependent receptor family" AND ' \
        'locations:(location:"Cell outer membrane [SL-0040]") ' \
        'NOT partial NOT fragment','format':'fasta'}
    print("Fetching TonB-depenent receptors")
    fetch(url, params, dbpath+'iromps.faa')
    return dbpath+'iromps.faa'

def run_blast(in_file, length, cpus, evalue, blastdb, fastadb):
     cmd = ['blastp', '-db', dbpath+blastdb,
           '-outfmt', '5', '-max_hsps', '1',
           '-evalue', evalue, '-max_target_seqs', '10000',
           '-num_threads', cpus, '-query', '-']
     blast = Popen(cmd, stdout=PIPE, stdin=PIPE)
     blast.stdin.write(in_file.encode())
     blast_out = StringIO(blast.communicate()[0].decode('utf-8'))
     hit_acc = []
     for q in SearchIO.parse(blast_out, "blast-xml"):
        if len(q.hits) > 0:
            for hsp in q.hsps:
                if length-3 <= hsp.hit_end <= length+3:
                    hit_acc.append(hsp.hit_id)
     hits = ''
     for r in SeqIO.parse(fastadb, 'fasta'):
        if r.id in hit_acc:
            hits = hits + r.format("fasta")
     return hits
             
def align(in_file):
    cmd = ['muscle']
    muscle = Popen(cmd, stdout=PIPE, stdin=PIPE)
    muscle.stdin.write(in_file.encode())
    return StringIO(muscle.communicate()[0].decode('utf-8'))
    
             

def make_db_df(in_file):
    print("Creating DataFrame...")
    data = []
    for r in list(SeqIO.parse(in_file, "fasta")):
        if re.search(r'GN=(.*) PE=', r.description):
            gene = re.search(r'GN=(.*) PE=', r.description).group(1)
        else:
            gene = 'NA'
        data.append(r.id+','+gene+','+
                    re.search(r' (.*) OS=', r.description).group(1).replace(",", " ")+','+
                    re.search(r'OS=(.*) OX=', r.description).group(1)+','+
                    str(len(r.seq)))
    return pd.DataFrame([sub.split(",") for sub in data],
                              columns=['accession','gene','description','organism','length'])

def fetch_seed(acc):
    Entrez.email = "tomdstanton@gmail.com"
    return Entrez.efetch(
        db='protein',id=acc, rettype="fasta",
        retmode="text").read().rstrip('\n')


def main(): 
    if not os.path.isfile(dbpath+'iromps.faa'):
        fastadb = get_fastadb()
    else:
        fastadb = dbpath+'iromps.faa'

    fastadb = cluster(fastadb, 0.9)
    make_blast_db(fastadb, 'prot', 'irompdb')
    seed_df = pd.read_csv(dbpath+'iromps.csv')
    iromp = fetch_seed(seed_df['NCBI_Accession'][0])
    length = sum(not c.isspace() for c in iromp.split(']')[1])
    IroN = run_blast(iromp, length, '4', '1e-130', 'irompdb', fastadb)
    alignment = align(IroN)
    