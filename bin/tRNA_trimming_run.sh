#!/bin/bash

# Feb 6 2019
# Matt Berg

#get sample name
DB=`basename $sample | cut -d "." -f1`

perl ./bin/tRNA_trimming.pl blast_analysis_output_cov10/tRNA_sequence_$i.fasta $i

done
