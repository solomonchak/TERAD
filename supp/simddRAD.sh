#!/bin/bash

file="$1"

# 1. linearize sequence (genome fasta file separates into lines, that cause problem)
sed -e 's/\(^>.*$\)/#\1#/' "$file" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' >linear.fasta

# 2. Get all the contig, each in one line. (get every line that doesn't match '>')
grep -v '>' linear.fasta > contigs.txt
rm linear.fasta

#####################################################################
# 3. Cut with restriction enzymes: SbfI-EcoRI
	grep -o 'GAATTC.*' contigs.txt |grep -o '.*CCTGCAGG' > SbEc_digested_pt1.txt
	grep -o 'CTTAAG.*' contigs.txt |grep -o '.*GGACGTCC' > SbEc_digested_pt2.txt
	grep -o 'GGACGTCC.*' contigs.txt |grep -o '.*GAATTC' > SbEc_digested_pt3.txt
	grep -o 'CCTGCAGG.*' contigs.txt |grep -o '.*CTTAAG' > SbEc_digested_pt4.txt
	cat SbEc_digested_pt1.txt SbEc_digested_pt2.txt SbEc_digested_pt3.txt SbEc_digested_pt4.txt> SbEc_digested.txt
	rm SbEc_digested_pt1.txt SbEc_digested_pt2.txt SbEc_digested_pt3.txt SbEc_digested_pt4.txt

# 4. Size selection (‘‘wide’’ automated size selection (300 +/- 36 bp )
awk 'length>264' SbEc_digested.txt |  awk 'length<336' > "${1%.*}_SbEc_ddRAD.fa"
rm  SbEc_digested.txt

#####################################################################
# 5. Cut with restriction enzymes: SphI-EcoRI 
	grep -o 'GAATTC.*' contigs.txt |grep -o '.*GCATGC' > SpEc_digested_pt1.txt
	grep -o 'CTTAAG.*' contigs.txt |grep -o '.*CGTACG' > SpEc_digested_pt2.txt
	grep -o 'CGTACG.*' contigs.txt |grep -o '.*GAATTC' > SpEc_digested_pt3.txt
	grep -o 'GCATGC.*' contigs.txt |grep -o '.*CTTAAG' > SpEc_digested_pt4.txt
	cat SpEc_digested_pt1.txt SpEc_digested_pt2.txt SpEc_digested_pt3.txt SpEc_digested_pt4.txt> SpEc_digested.txt
	rm SpEc_digested_pt1.txt SpEc_digested_pt2.txt SpEc_digested_pt3.txt SpEc_digested_pt4.txt

# 6. Size selection (‘‘wide’’ automated size selection (300 +/- 36 bp )
awk 'length>264' SpEc_digested.txt |  awk 'length<336' > "${1%.*}_SpEc_ddRAD.fa"
rm  SpEc_digested.txt

#####################################################################
# 7. Cut with restriction enzymes:  EcoRI-MspI
	grep -o 'GAATTC.*' contigs.txt |grep -o '.*CCGG' > EM_digested_pt1.txt
	grep -o 'CTTAAG.*' contigs.txt |grep -o '.*GGCC' > EM_digested_pt2.txt
	grep -o 'GGCC.*' contigs.txt |grep -o '.*GAATTC' > EM_digested_pt3.txt
	grep -o 'CCGG.*' contigs.txt |grep -o '.*CTTAAG' > EM_digested_pt4.txt
	cat EM_digested_pt1.txt EM_digested_pt2.txt EM_digested_pt3.txt EM_digested_pt4.txt > EM_digested.txt
	rm EM_digested_pt1.txt EM_digested_pt2.txt EM_digested_pt3.txt EM_digested_pt4.txt

# 8. Size selection (‘‘wide’’ automated size selection (300 +/- 36 bp )
awk 'length>264' EM_digested.txt |  awk 'length<336' > "${1%.*}_EM_ddRAD.fa"
rm  EM_digested.txt

#####################################################################
# 9. Cut with restriction enzymes:  SphI-MluCI 
	grep -o 'AATT.*' contigs.txt |grep -o '.*GCATGC' > SM_digested_pt1.txt
	grep -o 'TTAA.*' contigs.txt |grep -o '.*CGTACG' > SM_digested_pt2.txt
	grep -o 'CGTACG.*' contigs.txt |grep -o '.*AATT' > SM_digested_pt3.txt
	grep -o 'GCATGC.*' contigs.txt |grep -o '.*TTAA' > SM_digested_pt4.txt
	cat SM_digested_pt1.txt SM_digested_pt2.txt SM_digested_pt3.txt SM_digested_pt4.txt > SM_digested.txt
	rm SM_digested_pt1.txt SM_digested_pt2.txt SM_digested_pt3.txt SM_digested_pt4.txt 

# 10. Size selection (‘‘wide’’ automated size selection (300 +/- 36 bp )
awk 'length>264' SM_digested.txt |  awk 'length<336' > "${1%.*}_SM_ddRAD.fa"
rm  SM_digested.txt

#####################################################################
# 11. Cut with restriction enzymes:    NlaIII-MluCI
	grep -o 'AATT.*' contigs.txt |grep -o '.*CATG' > NM_digested_pt1.txt
	grep -o 'TTAA.*' contigs.txt |grep -o '.*GTAC' > NM_digested_pt2.txt
	grep -o 'GTAC.*' contigs.txt |grep -o '.*AATT' > NM_digested_pt3.txt
	grep -o 'CATG.*' contigs.txt |grep -o '.*TTAA' > NM_digested_pt4.txt
	cat NM_digested_pt1.txt NM_digested_pt2.txt NM_digested_pt3.txt NM_digested_pt4.txt > NM_digested.txt
	rm NM_digested_pt1.txt NM_digested_pt2.txt NM_digested_pt3.txt NM_digested_pt4.txt 

# 12. Size selection (‘‘wide’’ automated size selection (300 +/- 36 bp )
awk 'length>264' NM_digested.txt |  awk 'length<336' > "${1%.*}_NM_ddRAD.fa"
rm  NM_digested.txt

rm contigs.txt
