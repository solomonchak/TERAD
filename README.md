# TERAD
Extraction of transposable element composition from RADseq data

## Installation

### 1. Download TERAD
Place contents in, for example, ~/Desktop/TERAD

###  2. Install RepeatMasker 
Install RepeatMasker (http://www.repeatmasker.org/RMDownload.html) and its dependable programs including RMBlast (http://www.repeatmasker.org/RMBlast.html) and TRF (http://tandem.bu.edu/trf/trf.html). 
Then run the RepeatMasker configuration script. 
Test run RepeatMasker and RepeatProteinMask before proceeding.
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
cd ~Desktop/TERAD
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

## Questions?
Email `solomonchak@gmail.com`
