#!/bin/bash

# tally_restriction_site.sh
# Use to tally RE sites in all enzymes and enzymes combination from Peterson et al. 2012 ddRAD paper.


file="$1"


# Total   
		
		echo Total; cat $file|grep '>'| wc -l

# MspI CCGG or GGCC 

	 echo MspI; cat $file| grep 'ccgg\|ggcc' | wc -l 
   
   
# EcoRI GAATTC or CTTAAG 

	echo EcoRI; cat $file| grep 'gaattc\|cttaag' | wc -l
   
   
# SphI GCATGC\|CGTACG

	   echo SphI; cat $file| grep 'gcatgc\|cgtacg' | wc -l
   
   
# MluCI AATT TTAA

	   echo MluCI; cat $file| grep 'aatt\|ttaa' | wc -l
      

# NlaIII CATG GTAC
   		
   		echo NlaIII; cat $file| grep 'catg\|gatc' | wc -l
   
   
# SbfI   CCTGCAGG GGACGTCC
   
   		echo SbfI; cat $file| grep 'cctgcagg\|ggacgtcc' | wc -l
   		
# HindIII	AAGCTT TTCGAA
   	
   		echo SbfI; cat $file| grep 'aagctt\|ttcgaa' | wc -l

# EcoRI-MspI

	echo EcoRI-MspI; cat $file| grep 'gaattc\|cttaag' | grep 'ccgg\|ggcc' | wc -l

# SbfI-EcoRI  

		echo SbfI-EcoRI; cat $file| grep 'cctgcagg\|ggacgtcc' | grep 'gaattc\|cttaag' | wc -l

# SphI-EcoRI     

		echo SphI-EcoRI; cat $file| grep 'gcatgc\|cgtacg' | grep 'gaattc\|cttaag' | wc -l

# EcoRI-MspI 

		echo EcoRI-MspI; cat $file| grep 'gaattc\|cttaag' | grep 'ccgg\|ggcc' | wc -l

# SphI-MluCI     

		echo MspI; cat $file| grep 'gcatgc\|cgtacg' | grep 'aatt\|ttaa' | wc -l


# NlaIII-MluCI

		echo NlaIII-MluCI; cat $file| grep  'catg\|gatc' | grep 'aatt\|ttaa' | wc -l
		
# SphI - HindIII		


		echo SphI-HindIII; cat $file| grep  'gcatgc\|cgtacg' | grep 'aagctt\|ttcgaa' | wc -l









# END