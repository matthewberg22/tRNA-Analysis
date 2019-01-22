#!/bin/bash

# Jan 22 2019
# Matt Berg

#get sample name
DB=`basename $sample | cut -d "." -f1`

perl ./bin/BLAST_analysis.pl blast_output/$DB.blast $DB 10

done

# move to output folder
mv BLAST_analysis_*.txt blast_analysis_output_cov10/
mv tRNA_sequence_*.fasta blast_analysis_output_cov10/
rm Allcoverage_*.txt
