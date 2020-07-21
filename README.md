# SideroScanner
###### A tool for annotating IROMPs in bacteria
_By Tom Stanton_ \
Schneiders Lab - University of Edinburgh

Issues/queries/advice?
[email me!](T.D.Stanton@sms.ed.ac.uk) \
(Disclaimer: I'm not a bioinformatician so there are bound to be some bugs...)

[![alt text][1.1]][1]
[![alt text][6.1]][6]

[1]: http://twitter.com/tomstantonmicro
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)
[6]: http://www.github.com/tomdstanton
[6.1]: http://i.imgur.com/0o48UoR.png (github icon with padding)

### Introduction
*What does SideroScanner do?*
* Annotates IROMPs with accuracy, sensitvity and speed!
* IROMP = **I**ron **R**egulated **O**uter **M**embrane **P**rotein

* Aims to provide as much biologically relevant information as possible for
downstream wet/dry-lab validation.

*Why do we care about these proteins?*
* They are important **virulence factors**.
* They are the entry point of
**siderophore-antibiotic conjugates** and
have a role in mediating **resistance**
to these compounds.

*How does it work?*
1. Takes AA or DNA fasta (translates DNA to AA).
2. Iteratively filters proteins with conserved IROMP domain architecture.
3. Annotates them with a curated HMM library.

*Extra features for assemblies:*
* Screens for **plasmids/MGEs** to determine genomic location of hits.
* Screens hit flanking genes for additional **virulence factors**.
* Screens hit promoter regions for **Fur binding sites**.

**Please cite:**
```
SideroScanner: A tool for annotating IROMPs in bacteria.
Thomas David Stanton, 2020
https://github.com/tomdstanton/sideroscanner
```
### Dependencies
All of the following can be installed with conda:

|   sideroscanner  | sideroscanner-buildhmms | sideroscanner-builddbs |
|:----------------:|:-----------------------:|:----------------------:|
|   python >=3.7   |       python >=3.7      |      python >=3.7      |
| prodigal >=2.6.3 |       hmmer >=3.3       |      blast+ >=2.9      |
|    hmmer >=3.3   |       trimal >=1.4      |     cd-hit >=4.8.1     |
|   orfm >=0.7.1   |       muscle >=3.8      |                        |
|   **optional:**  |       blast+ >=2.9      |                        |
|   blast+ >=2.9   |      cd-hit >=4.8.1     |                        |
|   meme >=5.0.5   |                         |                        |

### Installation
Clone repo and install with python:
```
git clone --recursive https://github.com/tomdstanton/sideroscanner
cd sideroscanner
python setup.py install
```
## Usage
### Annotate
**Example commands:**
* Scan some gzipped fasta files and export results: \
```sideroscanner -i path/to/*.fna.gz -o results.csv -e annotated_sideroscanner.faa```
* Annotate a protein from NCBI: \
```efetch -db protein -id "WP_004151913.1" -format fasta | sideroscanner -i -```
* Scan an assembly, determine flanking genes, genomic location of hits
and putative Fur-binding-sites 300bp upstream of hits: \
```sideroscanner -i genome.fna -b 300 -f -l```

**Example output:**

| query          | hit   | description   |   score |   len |   kDa |
|:---------------|:------|:--------------|--------:|------:|------:|
| WP_004151913.1 | Fiu   | DHBS receptor |    1314 |   761 |  81.6 |

* Results are printed in markdown format so they can be pasted easily
into a variety of different programs and it looks nice(ish)! However
this can be changed if there is a demand for more pipe-able features.
* Takes any number of protein or DNA fasta-formatted files, accepts gzipped files too.
The program iterates over input files but not individual records in the file, 
for example if inputting a file of contigs, it will assume the whole file 
is a *sample* as opposed to treating each record as an individual *sample*.
```
Options:
  -i [- [- ...]]  path/to/input ['-' for STDIN]
                  -----------------------------------------------
  -o [-]          output results.csv instead of STDOUT
                  [optional: path/to/output/file]
                  [default: sideroscanner_DDMMYY_hhmm.csv]
                  -----------------------------------------------
  -l [int]        genomic location screen
                  [optional: blastn percid] [default: 90]
                  -----------------------------------------------
  -f [int]        flanking CDS screen
                  [optional: number of up/downstream CDS]
                  [default: 3]
                  -----------------------------------------------
  -b [int]        Fur binding site screen
                  [optional: length of promoter region]
                  [default: 200]
                  -----------------------------------------------
  -e [-]          export annotated proteins
                  [optional: path/to/export/fasta]
                  [default: filename_sideroscanner.faa]
                  -----------------------------------------------
  -t int          number of threads [default: max cpus]
                  -----------------------------------------------
  --lowqual [-]   'meta' CDS prediction AND single domain filter
                  [optional: path/to/draft/genome]
                  [default: all inputs]
                  -----------------------------------------------
  --nofilter [-]  turn off domain filter
                  [optional: path/to/input/fasta]
                  [default: all inputs]
                  -----------------------------------------------
  --lib hmm       path/to/custom/HMM/library.hmm
                  -----------------------------------------------
  -v              print version and exit
                  -----------------------------------------------
  --logo          super secret surprise!
                  -----------------------------------------------
  -h              show this help message and exit
```
### Build databases
By default, SideroScanner comes with just the IROMP HMM library.
* The hit location command (```-l```) uses [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb/)
and [ICEBerg2.0](https://db-mml.sjtu.edu.cn/ICEberg/),
however it will work if either one is absent: \
```sideroscanner-builddbs -db plsdb mgedb```
* The flank command (```-f```) uses a concatenated, non-redundant database from
[Patric-VF](https://www.patricbrc.org/), [Victors](http://www.phidias.us/victors/index.php)
and [VFDB](http://www.mgc.ac.cn/VFs/main.htm): \
```sideroscanner-builddbs -db flankdb```
* The Fur binding site command (```-b```) uses a concatenated 
Fur-box pwm from [CollecTF](http://www.collectf.org/browse/home/) which
unfortunately cannot be downloaded, but comes with SideroScanner.
```
Options:
  -db plsdb/mgedb/flankdb [plsdb/mgedb/flankdb ...]
                        build specific database [default: all]
                        -----------------------------------------------
  --keep                keep database fasta files
                        -----------------------------------------------
  -h                    show this help message and exit
```
* There might be some bugs if either the plsdb/mgedb is downloaded without the other
so it safer to download both.
### Build HMM library
* Adding new HMMs to the library is easy! Just run the following command
and follow the promts to enter a protein name, **valid** NCBI accession
and a short description:\
```sideroscanner-buildhmms --append```
* You then have the option of either appending the new HMM to the library or 
rebuilding it from scratch with your chosen options.
```
Options:
  -o [-]    save info about proteins in seed alignment
            [optional: path/to/output/file]
            [default: seed_alignment_DDMMYY_hhmm.xlsx]
            -----------------------------------------------
  -w int    protein length window [default: 2]
            -----------------------------------------------
  -e str    evalue for blastp [default: 1e-130]
            -----------------------------------------------
  -t int    number of threads [default: max cpus]
            -----------------------------------------------
  --keep    keep seed alignments
            -----------------------------------------------
  --review  print table of reference proteins for inspection
            -----------------------------------------------
  --append  add new protein to reference table
            -----------------------------------------------
  -h        show this help message and exit
```
![Image](https://github.com/tomdstanton/sideroscanner/blob/master/sideroscanner.png)