#!/bin/bash

# By Matt Berg
# May 2019

#Script runs extract_analysis.pl for all individual samples in the dataset

echo -e "Group	Identifier	TotaltRNAs	TotalReads	ConfidentReads	AmbiguousReads	DuplicatedtRNAs" > ./data/summarized_output_Cov10.txt

for i in GL4113 GL4181 GL16677 GL4614 GL195 GL4296 GL4290 GL4265 GL4264 GL4263 GL4262 GL4261 GL1481 GL996 GL594 GL16678 GL2100 GL1609 GL16679 GL15191 GL12471 GL4865 GL4657 GL150 GL4277 GL4253 GL2157 GL151 GL141 GL525 GL503 GL497 GL481 GL480 GL475 GL378 GL336 GL244 GL499 GL328 GL192 GL597 GL452 GL191 GL4297 GL16711 GL16714 GL16715
do
echo -n "Control	" >> ./data/summarized_output_Cov10.txt
perl ./bin/extract_analysis.pl ./data/blast_analysis_output_cov10/BLAST_analysis_$i.txt >> ./data/summarized_output_Cov10.txt
done

for i in GL4937 GL4938 GL6311 GL6453 GL6477 GL9006 GL11601 GL11620 GL11962 GL12586 GL12231 GL12932 GL3346 GL3755 GL3760 GL3952 GL4047 GL4147 GL4667 GL4869 GL4912 GL9076 GL12511 GL12758 GL13 GL157 GL193 GL434 GL485 GL585 GL589 GL699 GL710 GL769 GL836 GL1026 GL1312 GL1327 GL1377 GL1511 GL1575 GL1590 GL1845 GL1918 GL2234 GL2604 GL6488 GL11637
do
echo -n "HTG	" >> ./data/summarized_output_Cov10.txt
perl ./bin/extract_analysis.pl ./data/blast_analysis_output_cov10/BLAST_analysis_$i.txt >> ./data/summarized_output_Cov10.txt
done
