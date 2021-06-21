# SideroScanner
###### A tool for annotating IROMPs in bacteria
_By Tom Stanton_ :scientist: \
[![alt text][1.1]][1] [![alt text][6.1]][6] \
Issues/queries/advice?
[email me!](mailto:s1895738@ed.ac.uk?subject=[sideroscanner])

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
1. Takes AA or DNA fasta file or input stream.
2. Iteratively filters proteins with conserved IROMP domain architecture.
3. Annotates them with a curated HMM library.

*Extra features for assemblies:*
* Screens for **plasmids/MGEs** to determine genomic location of hits.
* Screens hit flanking genes for additional **virulence factors**.
* Screens hit promoter regions for **Fur binding sites**.
* Parameter adjustment for low quality assemblies.

**Please cite:**
```
SideroScanner: A tool for annotating IROMPs in bacteria.
Thomas David Stanton, 2020
https://github.com/tomdstanton/sideroscanner
```
### Dependencies
```
python >=3.8
meme >=5.0.5 (older meme versions have a different output and will break!)
prodigal >=2.6.3
hmmer >=3.3
blast ~=2.9.0
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
```sideroscanner *.fna.gz -o results.csv -e annotated_sideroscanner.faa```
* Find IROMPs in an assembly, their genomic location and flanking genes: \
```sideroscanner assembly.fna -l -f```
* Annotate a protein from NCBI: \
```efetch -db protein -id "WP_004151913.1" -format fasta | sideroscanner```
* Find IROMPs in an NCBI genome assembly and predict Fur binding sites
500bp upstream: \
```wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz | gunzip -c | sideroscanner -b 500 ```

```
usage: sideroscanner <IN.fasta> [options]
Options:
  -o [-]          output results.csv instead of markdown STDOUT
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
  -t int          number of threads [default: all]
                  -----------------------------------------------
  --lowqual [-]   'meta' CDS prediction AND single domain filter
                  [optional: path/to/draft/genome]
                  [default: all inputs]
                  -----------------------------------------------
  --nofilter [-]  turn off domain filter
                  [optional: path/to/input/fasta]
                  [default: all inputs]
                  -----------------------------------------------
  --molecule      force molecule type: genome / gene / protein
                  [default: guess molecule type]
                  -----------------------------------------------
  --lib hmm       path/to/custom/HMM/library.hmm
                  -----------------------------------------------
  -v              print version and exit
                  -----------------------------------------------
  --logo          super secret surprise!
                  -----------------------------------------------
  -h              show this help message and exit
```
**Known issues/bugs:**
* Running makeblastdb in blast v2.10 using ```sideroscanner-builddbs```
  can crash because of the insane amount of memory it tries to use.
* SideroScanner can get confused between the 
  highly similar enterobactin/salmochelin receptors FepA, PfeA and IroN.
  It's worth checking out the flanking proteins with ```-f``` to determine
  what  they really are.
* If the flanking proteins are too long, this can cause blastp to hang.
This has been observed in *P.aeruginosa* PAO1.
* Gzipped files cannot be piped in because I haven't figured out how to
decompress STDIN in python. Just use ```gunzip -c | sideroscanner```

**Features in development:**
* Outputting/appending GFF files.
* Further optimising HMM construction.
* Gzipped stdin.
### Build databases
By default, SideroScanner comes with just the IROMP HMM library.
* The hit location command (```-l```) uses [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb/)
and [ICEBerg2.0](https://db-mml.sjtu.edu.cn/ICEberg/),
however it will work if either one is absent: \
```sideroscanner-builddbs -db plsdb mgedb```
* The flank command (```-f```) uses a concatenated, non-redundant database from
[Patric-VF](https://www.patricbrc.org/), [Victors](http://www.phidias.us/victors/index.php),
[VFDB](http://www.mgc.ac.cn/VFs/main.htm) and 'siderophore'-related 
proteins from [UniProt](https://www.uniprot.org/uniprot/?query=siderophore+AND+taxonomy%3A%22Bacteria+%5B2%5D%22+NOT+receptor+NOT+partial+NOT+fragment&sort=score): \
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
  --overwrite           overwrite pre-existing databases
                        -----------------------------------------------
  -h                    show this help message and exit
```
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
  --report [-]  save info about proteins in seed alignment
                [optional: path/to/output/file]
                [default: seed_alignment_DDMMYY_hhmm.xlsx]
                -----------------------------------------------
  -w int        protein length window [default: 3]
                -----------------------------------------------
  -e str        evalue for blastp [default: 1e-130]
                -----------------------------------------------
  -t int        number of threads [default: 4]
                -----------------------------------------------
  --keep        keep seed alignments
                -----------------------------------------------
  --inspect     print table of reference proteins
                -----------------------------------------------
  --append      add new protein to reference table
                -----------------------------------------------
  -h            show this help message and exit
```
![Image](https://github.com/tomdstanton/sideroscanner/blob/master/sideroscanner.png)
