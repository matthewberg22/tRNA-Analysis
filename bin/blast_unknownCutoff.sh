#!/bin/bash

# May 2019
# Daniel Giguere

#blast trnas against each database
#set directory to top of project file
# usage: nohup ./blast.sh &

for sample in ./data/blast_output/*.blast; do

#get sample name
DB=`basename $sample | cut -d "." -f1`

# max target seqs set to 10 000 to ensure all reads are identified
blastn -query GtRNAdb_20flanking.fasta -db all_blast_dbs/$DB -perc_identity 95 -outfmt "7 std btop sseq qlen" -max_target_seqs 10000 -num_threads 40 > blast_output/$DB.blast

# arg 0 = blast output
# arg 1 = patient identifier
# arg 2 = coverage cutoff (0 if you don't know).
perl ./bin/BLAST_analysis.pl ./data/blast_output/$DB.blast $DB 0

done

# move to output folder
mv BLAST_analysis_*.txt ./data/blast_analysis_output/
mv tRNA_sequence*.fasta ./data/blast_analysis_output/
cat *Allcoverage_*.txt >> TotalCoverage.txt
rm Allcoverage_*.txt
mv TotalCoverage.txt ./data/blast_analysis_output/
