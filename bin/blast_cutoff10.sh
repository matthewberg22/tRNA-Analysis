#!/bin/bash

# May 2019
# Matt Berg

for sample in ./data/blast_output/*.blast; do

#get sample name
DB=`basename $sample | cut -d "." -f1`

perl ./bin/BLAST_analysis.pl ./data/blast_output/$DB.blast $DB 10

done

# move to output folder
mkdir ./data/blast_analysis_output_cov10
mv BLAST_analysis_*.txt ./data/blast_analysis_output_cov10/
mv tRNA_sequence*.fasta ./data/blast_analysis_output_cov10/
rm Allcoverage_*.txt
