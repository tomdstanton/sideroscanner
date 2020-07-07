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
### Annotate
Takes any number of protein or DNA fasta-formatted files, accepts gzipped files too.
The program iterates over input files but not individual records in the file, 
for example if inputting a file of contigs, it will assume the whole file 
is a sample as opposed to treating each record as an individual sample.
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

**Example commands:**
* Scan some gzipped fasta files and export results: ```sideroscanner.py -i path/to/*.fna.gz -o results.csv```
* Scan some fasta files with a bad egg: ```sideroscanner.py -i path/to/*.fasta --lowqual path/to/bad_egg.fasta```
* Scan an assembly and determine flanking genes plus genomic location of hits: ```sideroscanner.py -i genome.fna -f -l```
* Scan an assembly and export annotated proteins: ```sideroscanner.py -i genome.fna -e hits.faa```

**Example output:**
18 seconds on a 4-core laptop
```

```
### Build databases
By default, SideroScanner comes with just the IROMP HMM library because not everyone
may want to screen for plasmids of virulence factors. But those that do can just
run ```builddbs.py``` and/or specify the databases for the analysis of interest.
The ```sideroscanner.py -l``` command to determine location of hits uses both the
```plsdb``` from [PLSDB][https://ccb-microbe.cs.uni-saarland.de/plsdb/],
and ```mgedb```  from [ICEBerg][https://db-mml.sjtu.edu.cn/ICEberg/],
however if either one is absent
the analysis can still be performed with the other, depending on what you want to look for.
```
builddbs.py

Options:
  -db [{plsdb,mgedb,flankdb}]
                        | choose a specific database to build
                            [default: all]
                        -----------------------------------------------
  -h                    | show this help message and exit
```
### Build HMM library
By


![Image](https://github.com/tomdstanton/sideroscanner/blob/master/sideroscanner.png)
