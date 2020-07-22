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
```
meme >=5.0.5
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
* Scan an assembly and look for putative Fur-binding-sites 500bp upstream of hits: \
```sideroscanner -i NTUH-K2044_genomic.fna -b 500```

**Example output of the above command:**

| query            | hit           | description             |   score |   len |   kDa |   start |     end |   str | fur_start   | fur_end   | p_value   | fur_box                                            |
|:-----------------|:--------------|:------------------------|--------:|------:|------:|--------:|--------:|------:|:------------|:----------|:----------|:---------------------------------------------------|
| NC_012731.1_3308 | CirA          | DHBS/colicin receptor   |  1263.2 |   657 |  73   | 3648572 | 3650545 |    -1 | 3650715     | 3650765   | 8.16e-05  | TCATAATTGTATTGATAATCGTAATCATTATCTTTATCATTTTCGCCCAT |
| NC_012731.1_2780 | FitA          | Coprogen receptor       |  1158.3 |   690 |  75.6 | 3033742 | 3035814 |     1 | 3033280     | 3033330   | 9.24e-05  | CTGCAGTCAAATCAAAATAAATCTCATTCTCTTTTGATCATGACGGGGAT |
| NC_012731.1_857  | FhuA          | Ferrichrome receptor    |  1083.7 |   735 |  81.4 |  981000 |  983207 |     1 | 980937      | 980987    | 2.41e-05  | ATCGCCCGTCATAATAATAATTCTCGTTTACGTTATCATTCACTTTCATC |
| NC_012731.1_2580 | YncD          | TonB-dependent receptor |  1171.9 |   701 |  77.1 | 2835886 | 2837991 |     1 | -           | -         | -         | -                                                  |
| NC_012731.1_1970 | PfeA_ortholog | Enterobactin receptor   |  1310   |   727 |  80.1 | 2188602 | 2190785 |    -1 | 2190821     | 2190871   | 1.08e-05  | GAGAATATTAATGATAACAATTATCATTACAATGTAACGAGATGAATCTC |
| NC_012731.1_1976 | Fiu           | DHBS receptor           |  1314.7 |   761 |  81.5 | 2194552 | 2196837 |    -1 | 2197257     | 2197307   | 9.83e-06  | TGGTGATGGCCCCGGGATGGATCCGCACCGAGCTCGGCGGGGCAGATGCC |
| NC_012731.1_1322 | FepA          | Enterobactin receptor   |  1264.5 |   742 |  82.3 | 1489016 | 1491244 |    -1 | 1491398     | 1491448   | 2.24e-06  | TGATAATATTATTGATAACTATTTGCATTTGCAATAGCGTATTAGCGCGC |
| NC_012731.1_269  | FepA          | Enterobactin receptor   |  1215.9 |   752 |  82.7 |  318496 |  320754 |     1 | 318189      | 318239    | 3.15e-05  | CGGTAATAAAAATGAGATTCATTATCAAGATGATAATAATCAATATCGGA |
| NC_012731.1_3133 | IroN          | Salmochelin receptor    |  1365   |   724 |  79   | 3433588 | 3435762 |    -1 | -           | -         | -         | -                                                  |
| NC_006625.1_15   | IroN          | Salmochelin receptor    |  1333.2 |   724 |  79.3 |   13666 |   15840 |    -1 | 15930       | 15980     | 1.31e-05  | ATCGTAAGTAATGATAATTATTATCATTTGTGGGGAAGAAATTCAACCCT |
| NC_012731.1_2971 | FoxA          | Ferrioxamine receptor   |  1291.2 |   706 |  77.4 | 3248320 | 3250440 |     1 | -           | -         | -         | -                                                  |
| NC_012731.1_3123 | FyuA          | Yersiniabactin receptor |  1488.2 |   673 |  73.7 | 3424261 | 3426282 |     1 | -           | -         | -         | -                                                  |
| NC_006625.1_25   | FecA          | Ferric-citrate receptor |  1487.3 |   708 |  78.4 |   27193 |   29319 |     1 | -           | -         | -         | -                                                  |
| NC_012731.1_106  | BtuB          | Vitamin B12 receptor    |  1055.6 |   612 |  68.1 |  117391 |  119229 |     1 | 117196      | 117246    | 8.51e-05  | GTCTGCACCGAGTGCGCTCATTTCTCACCTTCCTTACCGCTGCGCGTCAG |
| NC_012731.1_3784 | HmuR          | Hemin receptor          |  1570.7 |   787 |  86.3 | 4179655 | 4182018 |     1 | 4179509     | 4179559   | 3.94e-05  | TTTGTTGCAAACGATAATACCTATCATTACCATTCGCAATCAACAAATGG |
| NC_012731.1_1250 | FcuA          | Ferrichrome receptor    |   811   |   732 |  79.2 | 1406993 | 1409191 |     1 | -           | -         | -         | -                                                  |
| NC_006625.1_214  | IutA          | Aerobactin receptor     |  1316.5 |   733 |  80.9 |  209679 |  211880 |    -1 | 212075      | 212125    | 7.51e-05  | GGGTTCGGCGATGCCATGGGTTTGCATACTGGCGTTGACCGCGAAGATGT |
| NC_012731.1_1819 | IutA          | Aerobactin receptor     |  1319.7 |   729 |  80.5 | 2043670 | 2045859 |    -1 | -           | -         | -         | -                                                  |

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