# Supporting scripts
For the paper TERAD: Extraction of transposable element composition from RADseq data (Chak et al. In revision)

## Supporting script 1: simddRAD.sh
In silico ddRADseq from genome assembly (simddRAD.sh)

This script will cut a genome assembly by the five sets of restriction enzymes from Peterson et al. 2012: SbfI-EcoRI, SphI-EcoRI, EcoRI-MspI, SphI-MluCI, and NlaIII-MluCI.

To run: 
`simddRAD.sh genome.fasta`

Outputs: 
- genome_SbEc_ddRAD.fa
- genome_SpEc_ddRAD.fa
- genome_EM_ddRAD.fa
- genome_SM_ddRAD.fa
- genome_NM_ddRAD.fa

## Supporting script 2: tallyRE.sh
In silico ddRADseq from genome assembly

Given a fasta file with multiple sequences (e.g., a TE database), this script will count the number sequences that will be cut by one of the restriction enzyme in Peterson et al. 2012 (EcoRI, MluCI, MspI, NlaIII, SbfI, and SphI), as well as the number of sequences that will be cut by two enzymes according to the enzyme combination in Peterson et al. 2012 (SbfI-EcoRI, SphI-EcoRI, EcoRI-MspI, SphI-MluCI, and NlaIII-MluCI.)

Example: 
`awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < TE_database.fa > TE_database_SingleLine.fa`

`sh tally_restriction_site.sh TE_database_SingleLine.fa | awk '{printf "%s%s",$0,NR%2?"\t":RS}' > TE_database_SingleLine_SingleLine_EZ_tally.txt`
