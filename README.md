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
* Annotates IROMPs with accuracy, sensitvity and speed!
* IROMP = **I**ron **R**egulated **O**uter **M**embrane **P**rotein

* Aims to provide as much biologically relevent information as possible for
downstream wet/dry-lab validation.

*Why do we care about these proteins?*
* They are important **virulence factors**.
* They are the entry point of
**siderophore-antibiotic conjugates** and
have a role in mediating **resistance**
to these compounds.

*How does it work?*
1. Takes AA or DNA fasta, preferrebly an assembly.
2. Translates DNA to AA.
3. Filters proteins with conserved IROMP domain architchture.
4. Annotates them with a custom profile hmm library.

Extra features:
* Screens assemblies for **plasmids/MGEs** to determine genomic location of the IROMPs.
* Screens flanking genes for additional **virulence factors**.

**Read the preprint here:**

### Dependencies
**Main program:** \
python >= 3.7 \
pandas \
biopython \
prodigal \
hmmer >= 3.3

**Additional analysis and dbs:** \
requests \
blast >=2.9

**Adding to and building the HMM library:** \
cd-hit \
muscle \
trimal

## Usage
**Main annotation program:**
```
sideroscanner.py -i query.fasta

Options:
  -i [- [- ...]]  | path/to/(i)nput/fasta
                  -----------------------------------------------
  -o [-]          | (o)utput file.csv instead of STDOUT
                      [optional: path/to/(o)utput/file]
                      [default: sideroscanner_DDMMYY_hhmm.csv]
                  -----------------------------------------------
  -l [int]        | determine genomic (l)ocation of hits
                      [optional: blastn percid]
                      [default: 90]
                  -----------------------------------------------
  -f [int]        | (f)lanking CDS screen
                      [optional: number of up/downstream CDS]
                      [default: 3]
                  -----------------------------------------------
  -e [-]          | (e)xport annotated proteins
                      [optional: path/to/export/fasta]
                      [default: path/to/'input'_sideroscanner.faa]
                  -----------------------------------------------
  -t int          | number of (t)hreads to use
                      [default: max available]
                  -----------------------------------------------
  --lowqual [-]   | (1) 'meta' CDS prediction AND (2) filters with plug domain only
                      [optional: path/to/(low)/(qual)ity/input/fasta]
                      [default: all inputs]
                  -----------------------------------------------
  --lib hmm       | path/to/custom/HMM/(lib)rary.hmm
                      [default: ~/sideroscanner/databases/irompdb/iromps.hmm]
                  -----------------------------------------------
  --dbpath path   | path/to/db/
                      [default: ~/sideroscanner/databases/]
                  -----------------------------------------------
  -v              | show version and exit
                  -----------------------------------------------
  -h              | show this help message and exit
```


### Examples
* Scan some fasta files: ```sideroscanner.py -i path/to/*.fasta -o results.csv```
* Scan some fasta files with a bad egg: ```sideroscanner.py -i path/to/*.fasta --lowqual path/to/bad_egg.fasta```
* Scan an assembly and determine flanking genes plus genomic location of hits: ```sideroscanner.py -i genome.fna -f -l```
* Scan an assembly and export annotated proteins: ```sideroscanner.py -i genome.fna -e hits.faa```

**Example output:**
18 seconds on a 4-core laptop
```
------------------------------------------------------------------------------------------------------------------------------------------------
SideroScanner: 0.0.1
Your system is Linux
2020-06-21-17:51:11
Using 4 threads...
------------------------------------------------------------------------------------------------------------------------------------------------
DNA input detected: SGH10
Extracting proteins...
5251 proteins extracted
Running hmmsearch...
Filtered 17 proteins with TonB-dependent Receptor Plug Domain
Running hmmsearch...
Filtered 17 proteins with TonB dependent receptor
Running hmmscan...
Annotated 17 proteins
Calculating length and molecular weight...
| query              | hit           | description                         |   score |   len |   kDa |   start |     end | str   |
|:-------------------|:--------------|:------------------------------------|--------:|------:|------:|--------:|--------:|:------|
| NZ_CP025080.1_1512 | CirA          | Catecholate/colicin receptor        |  1189.2 |   589 |  65.7 | 1603768 | 1605537 | +     |
| NZ_CP025080.1_2159 | FitA          | TonB-dependent siderophore receptor |  1353.6 |   690 |  75.6 | 2350967 | 2353039 | -     |
| NZ_CP025080.1_4115 | FhuA          | Ferrichrome receptor                |  1211.5 |   735 |  81.4 | 4450561 | 4452768 | -     |
| NZ_CP025080.1_2360 | YncD          | TonB-dependent siderophore receptor |  1443.2 |   701 |  77.1 | 2547544 | 2549649 | -     |
| NZ_CP025080.1_2968 | PfeA_ortholog | Enterobactin receptor               |  1378.8 |   727 |  80.1 | 3191333 | 3193516 | +     |
| NZ_CP025080.1_2962 | Fiu           | Catecholate receptor                |  1558.5 |   761 |  81.5 | 3185281 | 3187566 | +     |
| NZ_CP025080.1_3649 | FepA          | Enterobactin receptor               |  1409.6 |   742 |  82.3 | 3941098 | 3943326 | +     |
| NZ_CP025080.1_4705 | FepA          | Enterobactin receptor               |  1334.2 |   752 |  82.7 | 5112835 | 5115093 | -     |
| NZ_CP025081.1_120  | IroN          | Salmochelin receptor                |  1453.7 |   724 |  79.3 |  119074 |  121248 | +     |
| NZ_CP025080.1_1898 | FoxA          | Ferrioxamine receptor               |  1415.1 |   706 |  77.4 | 2088327 | 2090447 | -     |
| NZ_CP025080.1_1697 | FyuA          | Yersiniabactin receptor             |  1489.1 |   673 |  73.7 | 1855865 | 1857886 | -     |
| NZ_CP025081.1_110  | FecA          | Ferric-citrate receptor             |  1487.3 |   708 |  78.4 |  105595 |  107721 | -     |
| NZ_CP025080.1_4881 | BtuB          | Vitamin B12 receptor                |   890.3 |   612 |  68.1 | 5331412 | 5333250 | -     |
| NZ_CP025080.1_984  | HmuR          | Hemin receptor                      |  1644.9 |   787 |  86.3 | 1031991 | 1034354 | -     |
| NZ_CP025080.1_3721 | FcuA          | Ferrichrome receptor                |  1522.6 |   732 |  79.2 | 4023151 | 4025349 | -     |
| NZ_CP025081.1_147  | IutA          | Aerobactin receptor                 |  1435.1 |   733 |  80.9 |  147114 |  149315 | +     |
| NZ_CP025080.1_3152 | IutA          | Aerobactin receptor                 |  1426.1 |   729 |  80.5 | 3387018 | 3389207 | +     |
```
![Image](https://github.com/tomdstanton/sideroscanner/blob/master/sideroscanner.png)
