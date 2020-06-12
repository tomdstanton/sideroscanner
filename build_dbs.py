import requests
from subprocess import run
from tqdm import tqdm
import os


def fetch(url, params, outfile):
    try:
        r = requests.get(url, params, stream=True)
        total_size = int(r.headers.get('content-length', 0))
        block_size = 1024
        t = tqdm(total=total_size, unit='iB', unit_scale=True, position=0, leave=True,)
        with open(outfile, 'wb') as f:
            for data in r.iter_content(block_size):
                t.update(len(data))
                f.write(data)
        t.close()
        if total_size != 0 and t.n != total_size:
            print("ERROR, something went wrong")
        r.raise_for_status()
    except requests.exceptions.HTTPError as err:
        raise SystemExit(err)
    return r.content


def cluster(in_file, perc):
    out_file = os.path.splitext(in_file)[0]+'_nr.faa'
    print('Clutsering ' + os.path.basename(in_file) + ' at '+str(perc*100)+'%...')
    try:
        run(['cd-hit', '-i', in_file, '-o', out_file, '-c', str(perc), '-T', '0'])
        os.remove(out_file + '.clstr')
        os.remove(in_file)
    except:
        raise SystemExit()
    return out_file

def make_blast_db(infile, molecule, name):
    cmd = ['makeblastdb', '-in', infile,
           '-title', name, '-parse_seqids',
           '-dbtype', molecule, '-out', name]
    run(cmd)
    



def get_mgedb():
    mge_file = dbpath + 'mgedb'
    with open(mge_file, mode='wb') as localfile:
        localfile.write(fetch('http://202.120.12.136/ICEberg2/download/ICE_aa_all.fas'), None)


def get_plsdb():
    plsdb_file = dbpath + 'plsdb'
    with open(plsdb_file, mode='wb') as localfile:
        localfile.write(fetch('https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip'), None)
    with zipfile.ZipFile(dbpath + 'plsdb', "r") as zip_ref:
        zip_ref.extractall(dbpath)
    os.remove(dbpath + 'plsdb')
    os.remove(dbpath + 'plsdb.tsv')
    os.remove(dbpath + 'plsdb_changes.tsv')
    os.remove(dbpath + 'plsdb.msh')
    os.remove(dbpath + 'plsdb.abr')
    os.remove(dbpath + 'plsdb.sim')
    os.remove(dbpath + 'README.md')

# def get_flankdb():
#     handle1 = fetch(url,).decode('utf-8').splitlines()
#     records = list(SeqIO.parse(handle1, "fasta"))
#
#     fasta_file = open(input_file, 'r')  # opens inputfile for reading
#     fasta_lines = fasta_file.read().split('>')[0:]  # splits each sequence by header
#     unique_sequences = open(output_file, 'w')  # opens outputfile for writing
#
#     def remove_complete_duplicates(fasta_lines):
#         outputlist = []  # creates an empty list uniquesequences
#         setofuniqsequence = set()  # assigns a set for sequences -set can only contain uniquesequences
#         for sequence in fasta_lines:  # for loop - for sequence list in sequence list:
#             if sequence not in setofuniqsequence:  # If sequence has not been encountered yet i.e unique:
#                 outputlist.append(sequence)  # ... add it to list.
#                 setofuniqsequence.add(sequence)  # ... add it to set.
#         return outputlist
#
#     result = remove_complete_duplicates(fasta_lines)  # Remove duplicates from the list and defines as 'result'
#     unique_sequences.write('>'.join(result))  # '>' replaces lost > symbols
#     unique_sequences.close()  # closes output file
