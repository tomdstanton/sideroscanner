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
conda install python requests pandas biopython blast diamond cd-hit prodigal hmmer
```
## Usage

Input has to be a nucleotide or amino acid fasta file and molecule type will be auto-detected.
Proteins are extracted from nucleotide inputs using Prodigal, an accurate bacterial gene-predction program.
Then the TonB dependent receptor HMM filters proteins belonging to this family and passes them to Diamond
which performs a fast local blastp alignment against the curated database. Additionally, nucleotide
inputs can be screened for plasmid contigs to determine which hits are plasmid-encoded.
If you don't want the additional runtime of gene-prediction, you choose to perform a trasnlated blastx
alignemnt, however this will result in substantially fewer hits.

The fastest method is to use CDS protein fasta inputs, which are then filtered and aligned as before,
but without the additional gene-prediction step. You can also choose to turn off the HMM pre-filter in
both methods. It is worth noting that CDS predicted with Prodigal from a nucleotide file results in
more hits than CDS predicted with GeneMark and PGAP, although this might not always be the case.
N.B. plasmid screening can only be performed on nucleotide input files.

### Example usage: prepare database
Prepare default database: ```sideroscanner.py --makedb```

Process your own protein fasta file into a compatible DB: ```sideroscanner.py --makedb -d proteins.faa```

### Example usage: scan
Scan nucleotide: ```sideroscanner.py -s genome.fna```

Scan nucleotide and determine genomic location of hits: ```sideroscanner.py -s genome.fna -g```

Scan nucleotide quckly but fewer hits: ```sideroscanner.py -s genome.fna -x```

Scan proteins: ```sideroscanner.py -s refseq_cds.fna```

Scan proteins and turn off HMM pre-filter: ```sideroscanner.py -s refseq_cds.fna -b```

Scan proteins with additional hit descriptions: ```sideroscanner.py -s refseq_cds.fna -a```

For plasmid screening, download the plsdb (https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/) 
place the blast db files in sideroscanner script directory.
(-plsdb.fna.nsq
-plsdb.fna.nhr
-plsdb.fna.nin)