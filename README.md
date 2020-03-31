# SideroScannner

## A tool for accuratley annotating siderophore uptake proteins in bacteria.
Tom Stanton - Schneiders Lab - University of Edinburgh
T.D.Stanton@sms.ed.ac.uk

## Introduction
The ability to sequester iron is one of the key virulence determinants for 
pathogentic bacteria. In the host, free iron is scarce so bacteria 
synthesise and secrete iron-chelating small molecules called siderophores.
Once they have quired their ferric payload, they are uptaken back into the
cell via ligand-gated transporters, powered by the TonB protein.
These TonB-Dependent Transporters (TBDTs) are highly homologous, and are 
frequently misannotated by automated annotation pipelines. This can result in
the fallacious inferrence of a pathogens virulence or drug susceptibility.
SideroScanner seeks to mitigate that problem by using a HMM prefilter and
aligning to curated protein database of TBDTs to annotate as accuratley,
quickly and automatically as possible.

## Install with Conda

### Set up environment
```
conda create -n sideroscanner
source activate sideroscanner
```
### Add channels
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
### Install dependencies
```
conda install python pandas biopython blast diamond cd-hit prodigal hmmer
```
### You're good to go!

usage: sideroscanner.py [-h] [--makedb] [-s [SEQ [SEQ ...]]] [-o OUT] [-d DB]
                        [-m HMM] [-x] [-b] [-a] [-g]

optional arguments:
  -h, --help            show this help message and exit
  --makedb              Setup TBDT DB and HMM profile
  -s [SEQ [SEQ ...]], --seq [SEQ [SEQ ...]]
                        Path to fasta input, autodetects DNA or Protein.
  -o OUT, --out OUT     Path to output (comma-separated), otherwise prints to STDOUT.
  -d DB, --db DB        Path to protein database:
                        If used with --makedb, will process your own protein fasta file into a compatible DB.
                        If used with scan, a diamond-formatted database is required.
  -m HMM, --hmm HMM     Path to HMM profile, must be pressed in hmmer format.
  -x, --blastx          Perform translated alignment (protein files will default to blastp).
                        WARNING: This is the fastest option for nucleotide files but will result in fewer hits.
  -b, --blast_only      Turns off HMM pre-filter.
  -a, --annot           Add descriptions to hits
  -g, --genloc          Turns on chromosome/plasmid detection: 
                        Download the plsdb (https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)
                        Place the blast db files in sideroscanner script directory.
                            -plsdb.fna.nsq
                            -plsdb.fna.nhr
                            -plsdb.fna.nin
