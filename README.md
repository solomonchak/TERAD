# TERAD: Extraction of transposable element composition from RADseq data

Transposable elements (TEs)—selfish DNA sequences that can move within the genome—comprise a large proportion of the genomes of many organisms. Although low-coverage whole genome sequencing can be used to survey TE composition, it is non-economical for species with large quantities of DNA. Here, we utilize restriction site associated DNA sequencing (RADSeq) as an alternative method to survey TE composition. 

In our paper (in revision), we demonstrate in silico that double digest restriction-site associated DNA sequencing (ddRADseq) markers contain the same TE compositions as whole genome assemblies across arthropods. Then, we show empirically using eight *Synalpheus* snapping shrimp species with large genomes that TE compositions from ddRADseq and low-coverage whole genome sequencing are comparable within and across species. 

This bioinformatic pipeline, TERAD, is used to extract TE compositions from RADseq data.

*NOTE:* Our pipeline used only one end of the pair-end reads to remove the bias from the rarity of EcoRI cut sites among known Arthropod TEs. We found that the cut frequency of EcoRI was lower than that of MspI (56% vs. 85%, respectively). Therefore, we analyzed only the EcoRI-ends of the paired-end reads to include only TEs that did not have an EcoRI restriction site 

## Installation

### 1. Download TERAD
Place contents in, for example, ~/Desktop/TERAD

###  2. Install RepeatMasker 
Install RepeatMasker (http://www.repeatmasker.org/RMDownload.html) and its dependable programs including:
- RMBlast (http://www.repeatmasker.org/RMBlast.html) 
- TRF (http://tandem.bu.edu/trf/trf.html). 

Then run the RepeatMasker configuration script. 

Our program has been tested using RepeatMasker,v 1.332 2017/04/17 with RMBlast 2.9.0-p1

### 3. Install cd-hit 
https://github.com/weizhongli/cdhit

Our program has been tested using CD-HIT version 4.7 (built on Jan 25 2018).

### 4. Install R and packages
https://www.r-project.org/ 

Once R is installed, in terminal start R by typing `R`, then in the R console, type the following:
`install.packages("readr"); install.packages("plyr"); install.packages("fitdistrplus")`

If you are installing on a HPC system, you may need to choose custom library location (for example see: https://www.osc.edu/resources/getting_started/howto/howto_install_local_r_packages)

### 5. Add paths for RepeatMasker and cd-hit to PATH.
#### To set path on, for example, a Mac
To edit the bash_profile: `nano $HOME/.bash_profile`
add these lines at the end depending on the location where you have installed these programs: 
```
export PATH=$PATH:/Users/solomon/Documents/programs/cdhit-master
export PATH=$PATH:/Users/solomon/Documents/programs/RepeatMasker
```
To save: type `Control+ O`, then `Y`, then `Control+ X`

#### To set path on, for example, HPC
`nano .bashrc`
Add these lines:
```
PATH=/XXXX/XXXX/RepeatMasker:$PATH
PATH=/XXXX/XXXX/cdhit:$PATH
```

### 6. Now TERAD should be ready to go.
Test run:
```
cd ~/Desktop/TERAD
./TERAD test_file.fasta 4 ./arthro_ES_ND_PV_classified.fa none
```

## How to run
Inputs are:
1. A file to search for TE in fasta format
2. The number of cores to use
3. The custom library to search for TE or "none"
4. Query species for RepeatMasker or "none"
** Use either inputs 3 or 4. 
 
Make sure your input file, TERAD, extract_cdhit2.R, and RAD_TE_summary.R in the same folder as well as your custom TE library in fasta format if you wish to use one.

Test run:
`./TERAD test_file.fasta 4 ./arthro_ES_ND_PV_classified.fa none`
 or 
`./TERAD test_file.fasta 4 none arthropods`

## Main output
The main output is `test_file.fasta.cd.RAD_TE.summary2`. 

It is a .csv file that summarizes the proportions of RAD tags for each major TE subclass. `.int` and `.ext` indicate tags whether the restriction enzyme cut sites (e.g.,EcoRI if we used EcoRI-ends of the paired-end ddRAD reads) were `internal` or `external` to the TE. We only analyzed TEs that *don't* have the EcoRI sites (i.e., the `.ext` proportions) to avoid the low cut frequency of EcoRI in known Arthropod TEs.

| Column names  | Notes |
| ------------- | ------------- |
| sample  | Sample name   |
| lib.depth  | Total number of RAD tag   |
| DNA  | DNA transposons   |
| RC  | Helitron   |
| LTR  | LTR retrotransposon   |
| LINE  | Long interspersed nuclear elements    |
| SINE  | Short interspersed elements   |
| RNA  | RNA   |
| UO  | Unknown/Other   |
| Sim  | Simple repeat (Microsatellites)  |
| Sat  | Satellites   |
| LC  | Low complexity repeats   |
| TE  | Total TE   |
| TE.int  | Total TE where restriction sites are internal to the TE (inside the TE)   |
| TE.ext  | Total TE where restriction sites are external to the TE (outside the TE)    |
| NoTE  | RAD tags that have no TEs   |

## Other output from different components of the pipeline
From cd-hit
- test_file.fasta.cd
- test_file.fasta.cd.clstr

From extract_cdhit2.R
- test_file.fasta.cd.summary2
- test_file.fasta.cd.summary1

From RepeatMasker
- test_file.fasta.cd.cat
- test_file.fasta.cd.out
- test_file.fasta.cd.tbl
- test_file.fasta.cd.masked
- test_file.fasta.cd.RM.progress

From RepeatProteinMask
- test_file.fasta.cd.masked.protM.progress
- test_file.fasta.cd.masked.annot
- test_file.fasta.cd.masked.masked
- test_file.fasta.cd.masked.rmsimple.cat.all

From RAD_TE_summary.R
- test_file.fasta.cd.RAD_TE.summary1
- test_file.fasta.cd.RAD_TE.summary2


## Questions?
Email `solomonchak@gmail.com`
