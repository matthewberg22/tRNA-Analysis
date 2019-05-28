#!/bin/bash

# May 2019
# Matt Berg

for sample in ./data/blast_output/*.blast; do


#get sample name
DB=`basename $sample | cut -d "." -f1`

perl ./bin/tRNA_trimming.pl ./data/blast_analysis_output_cov10/tRNA_sequence-$DB\-Cov10.fasta $DB

done

mv total_alleles.txt ./data/
mv tRNAmutants.txt ./data/
