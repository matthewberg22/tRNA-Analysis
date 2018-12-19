# Workflow for processing tRNA reads
Daniel Giguere Oct 2 2018

**Summary:**

Reads were checked for quality with fastqc. The first 10 and last 2 bases were originally trimmed due to unexpected variations in per base sequence content. This interfered with read overlapping so raw reads were used.

```
# assign variables
FASTQC='/Volumes/data/bin/FastQC/fastqc'
READS='/Volumes/data/trna/reads/'

# output in fastqc_output
$FASTQC $READS -o .
```

################################################################################

# Merge overlapped read pairs with Pandaseq

Merge overlapping raw reads. Pandaseq (pandaseq 2.11 <andre@masella.name>) was run with default settings on macOS. Distrubition of length overlap by sample is plotted with proportion of reads merged by sample.

```
################################################################################
# run pandaseq

# this needs to be done with UNTRIMMED READS
nohup ./bin/pandaseq.sh &

# since i added "FINISHED" after each pandaseq, the stats can be grepped out
grep -B 9 "FINISHED" nohup.out > pandaseq_stats.out

################################################################################

# make a dataframe of the total number of reads and the number of reads that passed merging for all samples.

# from pandaseq_stats.out, match third column for READS, get number, paste together with tab separation for all reads.
# same for OK (reads that passed merging).
# get name of sample from the FINISHED identifier I put in by pasting substring starting at the 10th character (i.e., removes "FINISHED.")
# print SAMPLE, samples names store in header
# print READS, number of reads stored in x
# print OK, number of reads that passed stored in y.
awk '$3 ~ "READS" {x = x$4"\t"} $3 ~ "OK" {y = y$4"\t"} $1 ~ "FINISHED" {name = substr($0, 10); header = header name"\t"} END {print "SAMPLE\t"header"\nREADS\t"x"\nOK\t"y}' pandaseq_stats.out > passed_merging.txt


################################################################################

# plot statistics for each sample

R

d <- read.table("passed_merging.txt", sep = "\t", row.names= 1, header = TRUE)
# remove last column that shouldn't really be there.
d$X <- NULL
d.t <- data.frame(t(d))

# calculate proportion of reads that passed merging
d.t$prop <- (d.t$OK / d.t$READS) * 100

# plot percentage of reads merged against sample ID
pdf("figs/percentage_reads_merged_all_samples.pdf")
plot(d.t$prop, ylim = c(0,100), main = "Proportion of read pairs merged", xlab = "Sample Number", ylab = "% merged")
dev.off()

# plot length of overlaps for each sample now
awk '$3 ~ "OVERLAPS" {$1=$2=$3=""; x = $0} $1 ~ "FINISHED" {name = substr($0, 10); print name x} END {exit}' pandaseq_stats.out > overlaps.txt

R

d <- read.table("overlaps.txt", header = FALSE)

#only the integer
e <- data.frame(t(d[,2:153]))
colnames(e) <- d$V1

# plot each distribution
pdf("figs/overlap_distribution_all_samples.pdf")
colours <- colorRampPalette(c("blue", "red"))(96)
plot(e[,1], type="n", ylim=c(0, 12000), ylab = "Number of read pairs", xlab = "Length of overlap")

for (i in seq(e)) {
    points(e[,i], col=colours[i], pch = 19, cex = 0.2)
}

legend("topright", col=c("blue", "red"), pch = 19, legend = c("Sample 1", "Sample 96"))

dev.off()
################################################################################
```

### Blast against db of reads

For each sample, a BLAST database of the merged reads needs to be made.

Before making blastdb, you need to remove the long string that is common to each of the sample. eg. M00388:388:000000000-BGJ2L from M00388:388:000000000-BGJ2L:1:1101:15326:1913:2. BLAST will cut the identifier to a limited number of characters.

```
for sample in ./*merged.fasta; do

# get the sample name (i.e., GL141 from GL141_001_R1...)
NAME=`basename $sample | cut -d "_" -f1`

sed 's/M00388:[0-9]*:000000000-[a-zA-Z0-9]*//' $sample > $NAME-merged-edited.fasta

done
```

Make the DBs for each file and BLAST them.

```
#start with just sample GL141
nohup makeblastdb -in GL141_S2_L001_R1_001_edited.merged.fasta -dbtype nucl -out all_blast_dbs/GL141 &

DB="all_blast_dbs/GL141"
nohup blastn -query 20181215_tRNASeq_Query -db $DB -perc_identity 95 -outfmt 7 -max_target_seqs 10000 -num_threads 40 > GL141.blast &

#the rest of the files
nohup ./bin/makeblastdbs.sh &
```
Then blast the tRNAs against each database.

```
nohup ./bin/blast.sh &
```
