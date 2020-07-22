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

**Example output of the above command:**

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
| PAO1   | NC_002516.2_3343 | OptR          | TonB-dependent receptor               | 1197.8 | 721 | 79.4  | 3655986 | 3658151 | -1  | -                                               | sensor histidine kinase/response regulator                                                        | 3658566   | 3658616 | 5.12E-05 | GCACGGGCCACCGCCAACGGCTCCTGCCCGGCCTTCGGCGCTTCTGGCAG |
| PAO1   | NC_002516.2_950  | PirA          | Catecholate receptor                  | 1150.4 | 742 | 80.9  | 1018230 | 1020458 | 1   | gacS                                            | GTP pyrophosphokinase                                                                             | 1018038   | 1018088 | 2.36E-06 | TTACTCTCTTTTGTTAATGATTATCATCCGTGCCGATCTCCTCGGCCGGC |
| PAO1   | NC_002516.2_2984 | OptE          | TonB-dependent receptor               | 895.1  | 718 | 80.3  | 3265848 | 3268004 | 1   | -                                               | -                                                                                                 | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_2114 | OptM          | TonB-dependent receptor               | 1486   | 880 | 95.3  | 2269542 | 2272184 | -1  | -                                               | -                                                                                                 | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_4257 | FptA          | Pyochelin receptor                    | 841.2  | 802 | 87.4  | 4663854 | 4666262 | 1   | -                                               | -                                                                                                 | 4663688   | 4663738 | 6.61E-05 | AAATGGAGTAAATGATAGAGATTATCGTTATCAGTAGTATTCGCAGGAAT |
| PAO1   | NC_002516.2_444  | OptJ          | TonB-dependent receptor               | 1330.7 | 730 | 81.3  | 484964  | 487156  | 1   | -                                               | -                                                                                                 | 484764    | 484814  | 4.12E-05 | CGCCCGCCACGCGCGAGCCGGTGCCGCCCCGGTACGTCTGGCCGCCGGCC |
| PAO1   | NC_002516.2_4310 | FptA          | Pyochelin receptor                    | 1297.8 | 720 | 80    | 4724639 | 4726801 | -1  | fptC fptX                                       | pchI pchH pchG                                                                                    | 4726820   | 4726870 | 3.41E-06 | GATATACGATTTGATAATGCTTATCATTATGAAGAAATCGTCTCCGGGGC |
| PAO1   | NC_002516.2_3987 | FecA          | Ferric-citrate receptor               | 1345.5 | 784 | 85.5  | 4368837 | 4371191 | 1   | -                                               | peptide chain release factor 3                                                                    | 4368449   | 4368499 | 8.51E-05 | TGGACTTCGACCGCGAGCGACTCGGCGCTCCGCAGCCGCTGCCGGTCGGC |
| PAO1   | NC_002516.2_3873 | OprC          | TonB-dependent copper receptor        | 1210.6 | 660 | 72.5  | 4247703 | 4249685 | -1  | -                                               | 2-isopropylmalate synthase                                                                        | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_3487 | HasR          | Heme receptor                         | 1514.3 | 891 | 97.9  | 3814661 | 3817336 | -1  | hasA hasD hasE                                  | hasS                                                                                              | 3817356   | 3817406 | 6.58E-07 | GTTTTGTGTAATAAGAATCAGTCTCGTTAAGAGAACGATGGAGGCAGGGA |
| PAO1   | NC_002516.2_1296 | BtuB          | Vitamin B12 receptor                  | 345.7  | 616 | 67.6  | 1381804 | 1383654 | 1   | -                                               | -                                                                                                 | 1381637   | 1381687 | 5.82E-05 | GAAGGCGCGGCTGGAAGCGTCCAGCGCTTCGCTCGCGAGCCCGGAGACCG |
| PAO1   | NC_002516.2_2342 | YncD          | TonB-dependent receptor               | 916.2  | 710 | 78.2  | 2516429 | 2518561 | -1  | -                                               | -                                                                                                 | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_2134 | OptF          | TonB-dependent receptor               | 1580.2 | 883 | 95.6  | 2298012 | 2300663 | 1   | -                                               | -                                                                                                 | 2297613   | 2297663 | 1.74E-05 | GCTGCTGCCGATGGGAACCCGCACCCGTGGCCGGGAGGCGTTGCGCGAAC |
| PAO1   | NC_002516.2_2451 | FpvA          | Pyoverdine receptor                   | 735.8  | 800 | 89.6  | 2655232 | 2657634 | 1   | pvdE pvdF pvdO                                  | pvdD pvdJ pvdI                                                                                    | 2655109   | 2655159 | 8.49E-06 | GTTGCGAATAAATGAGATCGCTTCTCAATTCTCTTAGTATGACGTTTTAC |
| PAO1   | NC_002516.2_4245 | FvbA          | Vibriobactin receptor                 | 1250.2 | 694 | 76.8  | 4650374 | 4652458 | -1  | -                                               | -                                                                                                 | 4652587   | 4652637 | 9.90E-08 | CCCACCGCCCCGGGGAATCGGTATCATGGGCGCCGCCCGCGCCACCGCAG |
| PAO1   | NC_002516.2_2101 | SppR          | Xenosiderophore receptor              | 2045.6 | 985 | 108.1 | 2250858 | 2253815 | 1   | -                                               | -                                                                                                 | 2250429   | 2250479 | 4.69E-05 | GTTCGAACGGCTGGAGGATGCGTACCAGTCGACGCCGGCGGATCGCCTCC |
| PAO1   | NC_002516.2_799  | ZnuD          | Zinc receptor                         | 1136.1 | 687 | 74.2  | 849256  | 851319  | -1  | ATP-dependent protease La                       | bifunctional proline dehydrogenase/pyrroline-5-carboxylate dehydrogenase sodium/proline symporter | 851626    | 851676  | 6.97E-06 | GTGGCGCGGGAATACGGATCAGTATGGGTCAGGCCCATTCCGCCGGCGCG |
| PAO1   | NC_002516.2_1390 | OptN          | Siderophore receptor                  | 1669.2 | 813 | 89.1  | 1476384 | 1478825 | 1   | -                                               | -                                                                                                 | 1476206   | 1476256 | 2.25E-07 | GCGACACGGCCACGGGAAGGATACTCGGCAGCAGGTCGAGCAGCCGGTCG |
| PAO1   | NC_002516.2_2646 | OptF          | TonB-dependent receptor               | 1564.8 | 884 | 95    | 2930808 | 2933462 | -1  | -                                               | -                                                                                                 | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_1327 | TdhA          | Heme receptor                         | 757    | 851 | 95.2  | 1411585 | 1414140 | 1   | -                                               | -                                                                                                 | 1411464   | 1411514 | 1.90E-06 | GTCAGCATTAATGAAAATATTTTTCATTCGCATATGTGGGTTTTTCTCCC |
| PAO1   | NC_002516.2_1951 | BauA          | Acinetobactin receptor                | 716.5  | 804 | 86.4  | 2081853 | 2084267 | -1  | -                                               | -                                                                                                 | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_5000 | HasR          | Heme receptor                         | 273.5  | 989 | 108.3 | 5491346 | 5494315 | 1   | -                                               | -                                                                                                 | 5491202   | 5491252 | 1.19E-05 | ACGCAAAATTAGTGATAATAGTTTTCATTTGCGATTGAATCCAGCTGACT |
| PAO1   | NC_002516.2_197  | OptP          | TonB-dependent receptor               | 1304.6 | 790 | 84.6  | 219172  | 221544  | 1   | -                                               | -                                                                                                 | -         | -       | -        | -                                                  |
| PAO1   | NC_002516.2_2389 | SftP          | Hexylsulfate receptor                 | 1352.8 | 789 | 85    | 2577150 | 2579519 | 1   | -                                               | -                                                                                                 | 2576885   | 2576935 | 1.74E-05 | CGCGGCGCTGGCGGTTAACCTTTCGCAATCGGCCTTCAGCCGCAGCATCC |
| PAO1   | NC_002516.2_4773 | IutA          | Aerobactin receptor                   | 1038.6 | 742 | 81    | 5243178 | 5245406 | 1   | GTP-dependent nucleic acid-binding protein EngD | -                                                                                                 | 5243099   | 5243149 | 2.77E-06 | TTCGAATCATTTGATAATCATTATCGATTTGTTTAGCTTTGCCGCCCATC |
| PAO1   | NC_002516.2_1643 | CirA_ortholog | BtuB-like receptor                    | 1267.6 | 702 | 79    | 1756489 | 1758597 | -1  | -                                               | -                                                                                                 | 1758946   | 1758996 | 8.51E-05 | GGGCAAGGCGACGAAGATCAGTTGGCAATCGGCCAGGGTGCGCTCCAGGT |

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