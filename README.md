# SideroScanner
###### A tool for annotating IROMPs in bacteria
_By Tom Stanton_ \
Schneiders Lab - University of Edinburgh

Issues/queries/advice?
[email me!](T.D.Stanton@sms.ed.ac.uk)

[![alt text][1.1]][1]
[![alt text][6.1]][6]

[1]: http://twitter.com/tomstantonmicro
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)
[6]: http://www.github.com/tomdstanton
[6.1]: http://i.imgur.com/0o48UoR.png (github icon with padding)

### Introduction
*What does SideroScanner do?*
* Annotates **I**ron **R**egulated **O**uter **M**embrane **P**roteins with accuracy, sensitvity and speed!
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
```
python ~=3.7.6
meme >=5.0.5 (older meme versions have a different output and will break!)
prodigal >=2.6.3
hmmer >=3.3
blast >=2.9.0
trimal >=1.4.1
orfm >=0.7.1
muscle ~=3.8.1551
cd-hit >=4.8.1
```
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
* Feelin' crazy? Scan an NCBI genome assembly, screen the flanking genes and 500 bp upstream for
Fur binding sites: \
```wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz | gunzip -c | sideroscanner -i - -f -b 500 ```

**Example output of the above command (truncated):**

| sample | query            | hit           | description                           | score  | len | kDa   | start   | end     | str | upstream                                        | downstream                                                                                        | fur_start | fur_end | p_value  | fur_box                                            |
|--------|------------------|---------------|---------------------------------------|--------|-----|-------|---------|---------|-----|-------------------------------------------------|---------------------------------------------------------------------------------------------------|-----------|---------|----------|----------------------------------------------------|
| PAO1   | NC_002516.2_4939 | CntO          | Pseudopaline receptor                 | 1046.3 | 708 | 79.1  | 5427716 | 5429842 | -1  | -                                               | arginine decarboxylase                                                                            | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_2517 | FoxA          | Ferrioxamine receptor                 | 918    | 820 | 90    | 2782764 | 2785226 | -1  | -                                               | -                                                                                                 | 2785292   | 2785342 | 8.08E-06 | GACAAAAATAATGAAAATAATTTTCAGTGCGTTTCGCTCCGGACGAAAGC |
| PAO1   | NC_002516.2_157  | PA0151        | TonB-dependent receptor               | 1429.8 | 795 | 86.7  | 171047  | 173434  | 1   | -                                               | -                                                                                                 | 170992    | 171042  | 8.28E-07 | GTGCTCATGACCGCGATTGATTCGCATTAGAAAAGCGACCGACACGAGCA |
| PAO1   | NC_002516.2_4809 | PhuR          | Heme/hemoglobin receptor              | 1572.2 | 764 | 84.7  | 5289217 | 5291511 | 1   | phuS phuT phuU                                  | -                                                                                                 | 5289111   | 5289161 | 1.44E-05 | CTGAATCCATTTGATAATTATTTGCATTAGCGTTTTTCTGGCAGTACCTT |
| PAO1   | NC_002516.2_481  | FhuA          | Ferrichrome receptor                  | 798.6  | 802 | 88.2  | 530029  | 532437  | -1  | -                                               | -                                                                                                 | 532442    | 532492  | 4.19E-08 | GGGCTCCTGGATGGAAACGAGTCTCAATGCCTCCTTGCCGGACGAGTCGG |
| PAO1   | NC_002516.2_1347 | PfuA          | Hydroxamate-type siderophore receptor | 1272.2 | 732 | 80.8  | 1433166 | 1435364 | 1   | -                                               | -                                                                                                 | 1433106   | 1433156 | 2.77E-06 | CATCTTGGTTATTGAGAATCATTGGCATTTGATTGATGGAGGGTTTTTTT |
| PAO1   | NC_002516.2_1964 | CirA          | DHBS/colicin receptor                 | 490    | 653 | 72.3  | 2097491 | 2099452 | 1   | -                                               | -                                                                                                 | 2097111   | 2097161 | 1.74E-05 | ACCCGCGCCGCCGGGAAGCGTTCGAGCAGGCGCCGGCCGAGGGAGCCCGG |
| PAO1   | NC_002516.2_4608 | PfuA          | Hydroxamate-type siderophore receptor | 794.9  | 753 | 82.3  | 5053616 | 5055877 | -1  | -                                               | ibeC                                                                                              | 5056024   | 5056074 | 4.42E-06 | AATGATTGCCAATGATATTGATTTGCATTGGACATGTAAAACCGCTAGAG |
| PAO1   | NC_002516.2_2746 | FepA          | Enterobactin receptor                 | 1232   | 746 | 81    | 3040242 | 3042482 | 1   | vgrG1b                                          | -                                                                                                 | 3040172   | 3040222 | 5.34E-05 | TTACTCTCAAATAACAATCAATATCATTTGTGATCTCTTGCATTTCGCTG |

* Results are printed in markdown format so they can be pasted easily
into a variety of different programs and it looks nice(ish)! However
this can be changed if there is a demand for more pipe-able features.
* Takes any number of protein or DNA fasta-formatted files, accepts gzipped files too, however it seems
they cannot be piped currently.
* The program iterates over input files but not individual records in the file, 
e.g. if inputting a file of contigs, it will assume the whole file 
is a *sample* as opposed to treating each record as an individual *sample*.
* Fastq files are supported, but this feature doesn't work well. Optimisation is planned in future.
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
and follow the prompts to enter a protein name, **valid** NCBI accession
and a short description:\
```sideroscanner-buildhmms --append```
* TonB-dependent receptors will be downloaded from UniProt (~350Mb) and
made into a non-redundant BLASTP db.
* New proteins are blasted against the db using a defined E-value and
size window.
* Hits are aligned, trimmed and built into pHMMs.
* You have the option of overwriting all HMMs or just appending your new proteins.
```
Options:
  --report [-]    save info about proteins in seed alignment
                  [optional: path/to/output/file]
                  [default: seed_alignment_DDMMYY_hhmm.xlsx]
                  -----------------------------------------------
  -w int          protein length window [default: 2]
                  -----------------------------------------------
  -e str          evalue for blastp [default: 1e-130]
                  -----------------------------------------------
  -t int          number of threads [default: max cpus]
                  -----------------------------------------------
  --keep          keep seed alignments
                  -----------------------------------------------
  --review        print table of reference proteins for inspection
                  -----------------------------------------------
  --append        add new protein to reference table
                  -----------------------------------------------
  -h              show this help message and exit
```
![Image](https://github.com/tomdstanton/sideroscanner/blob/master/sideroscanner.png)