#!/bin/bash

# Dec 19 2018
# Daniel Giguere

#blast trnas against each database
#set directory to top of project file
# usage: nohup ./blast.sh &

for sample in all_blast_dbs/*.nhr; do

#get sample name of database
DB=`basename $sample | cut -d "." -f1`

#test
echo $NAME

# max target seqs set to 10 000 to ensure all reads are identified
blastn -query 20181215_tRNASeq_Query.fasta -db all_blast_dbs/$DB -perc_identity 95 -outfmt 7 -max_target_seqs 10000 -num_threads 40 > blast_output/$DB.blast

done
