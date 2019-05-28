#!/bin/bash

# May 2019
# Daniel Giguere

# make blast db for each sample
# set directory to top of project file

# make directory for output if it doesn't exist
mkdir -p all_blast_dbs

for sample in merged_reads/*merged-edited.fasta; do

# get only sample name (i.e., GL141)
NAME=`basename $sample | cut -d "-" -f1`

#test
#echo $NAME

makeblastdb -in $sample -dbtype nucl -out all_blast_dbs/$NAME;

done
