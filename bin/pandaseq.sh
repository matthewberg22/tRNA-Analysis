#!/bin/bash

# merge read pairs with pandaseq
# nov 20 2018
# daniel giguere

# the stats will be outputted to the nohup.out file (10+ GB in size)
# usage: nohup ./pandaseq.sh &

# set variable to where reads are located
READS='/raw_reads_all'

# make output directory if it doesn't exist
mkdir -p pandaseq

for f in $READS/*R1*.fastq.gz; do

# get sample names from file
B=`basename $f`
NAME=`echo $B | cut -d "." -f1`

R='R2'
# replace R1 with R2 to specify reverse read file
REV=`echo ${B/R1/$R}`

#test 
#echo $f
#echo $REV

pandaseq -f $f -r $READS/$REV -w pandaseq/$NAME.merged.fasta

#add unique identifier to grep out stats for each sample later
DONE="FINISHED"
echo $DONE.$NAME

done
