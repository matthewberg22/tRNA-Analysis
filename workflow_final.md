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

# Merge overlapped read pairs with Pandaseq

Merge overlapping raw reads. Pandaseq (pandaseq 2.11 <andre@masella.name>) was run with default settings on macOS. Distrubition of length overlap by sample is plotted with proportion of reads merged by sample.

```
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

#Librarys to make the graphs and clean the data
library(ggplot2)
library(tidyr)
library(cowplot)
library(ggrepel)

#A list of sample groups
Groups <- read.table(file.choose(), sep = "\t", header = TRUE)

#Shortens the identifiers to just the GL number
d.t$Sample <- rownames(d.t)
Split <- d.t %>% separate(Sample, c("Sample", NA, NA, NA, NA), "_")
d.t$Sample <- Split$Sample

#Adds which group corresponds to which GL number
d.t$group <- ordered(d.t$Sample, levels = Groups$Sample, labels = Groups$Group)
d.t <- d.t[order(d.t$group),]

#Graphs the percent merged
PercentMerged <- ggplot(d.t, aes(order(d.t$group), prop, col = group, label = Sample)) + geom_point() + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + xlab('') + ylab('Percent Merged') + geom_text_repel(aes(label=ifelse(prop<80,as.character(Sample),'')), force =50) + coord_cartesian(ylim = c(30,100))
svg("figs/percentage_reads_merged_all_samples.svg")
PercentMerged
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
```

### Blast against db of reads

For each sample, a BLAST database of the merged reads needs to be made.

Before making blastdb, you need to remove the long string that is common to each of the sample. eg. M00388:388:000000000-BGJ2L from M00388:388:000000000-BGJ2L:1:1101:15326:1913:2. BLAST will cut the identifier to a limited number of characters.

```
for sample in ./*merged.fasta; do

# get the sample name (i.e., GL141 from GL141_001_R1...)
NAME=`basename $sample | cut -d "_" -f1`

# the regex is required since difference samples have different values
sed 's/M00388:[0-9]*:000000000-[a-zA-Z0-9]*//' $sample > $NAME-merged-edited.fasta

done
```

Make the DBs for each file and BLAST them.

```
#start with just sample GL141
nohup makeblastdb -in GL141_S2_L001_R1_001_edited.merged.fasta -dbtype nucl -out all_blast_dbs/GL141 &

DB="all_blast_dbs/GL141"
nohup blastn -query all20_tRNAquery.fasta -db $DB  -perc_identity 95 -outfmt "7 std btop sseq qlen"   -max_target_seqs 10000 -num_threads 1 > $DB.blast &

#the rest of the files
nohup ./bin/makeblastdbs.sh &
```
Then blast the tRNAs against each database and run the blast anaylsis script. 

```
nohup ./bin/blast.sh &
```

Next, determine the ideal cutoff by graphing total coverage for every unique tRNA sequence observed from all individuals sequenced.

```
R
#Packages needed to make the graphs
library(ggplot2)
library(cowplot)
library(compositions)

#Reads in coverage summary file that is the output from my perl script
all.coverage <- read.table(file.choose(), header = TRUE, sep = "\t")

#Plots the frequency of coverage for all unique sequences seen in the 96 individuals sequenced
#Limits the x-axis to 300 to scale plot better
#Red line indicates the mean coverage
#Blue line is my coverage cut off right now
p0 <- ggplot(all.coverage, aes(x = Coverage)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim=c(0, 300), ylim=c(0,1000)) + xlab("Coverage of tRNA allele") + 
  ylab("Frequency") + 
  geom_vline(aes(xintercept = 5), color = 'orange') + 
  geom_vline(aes(xintercept = 10), color = 'red') + 
  geom_vline(aes(xintercept = 15), color = 'blue') +
  theme(legend.position="none")

#Centered log ratio transformation of the data
clr <- all.coverage[which(all.coverage$Coverage == 5), 5]
clr5 <- clr[1]

clr <- all.coverage[which(all.coverage$Coverage == 10), 5]
clr10 <- clr[1]

clr <- all.coverage[which(all.coverage$Coverage == 15), 5]
clr15 <- clr[1]

all.coverage$CLR <- as.vector(clr((all.coverage$Coverage)))

p1 <- ggplot(all.coverage, aes(x = CLR)) + 
  geom_histogram(binwidth = 0.5) + 
  coord_cartesian(ylim=c(0,10000)) + xlab("clr(Coverage of tRNA allele)") + 
  ylab("Frequency") + 
  theme(legend.position="none") +
  geom_vline(aes_(xintercept = clr15), color = 'blue') +
  geom_vline(aes_(xintercept = clr10), color = 'red') +
  geom_vline(aes_(xintercept = clr5), color = 'orange')

plot_grid(p0, p1, labels = c("A", "B"))

svg("figs/CoveragePlots.svg")
CoveragePlots
dev.off()

```
After looking at the cutoff graphs, cutoff was selected as 10x coverage. Ran perl script again to get list of tRNAs with >10x coverage.

```
nohup ./bin/tRNA_sequences_cutoff10.sh &
```
Ran tRNA trimming perl script on all samples to get tRNA mutant sequences without 20 bp 5' flanking and to get total counts for each loci.

Ran script to calculate allele frequency of each mutant. Output is a non-redundant list of mutant tRNA sequences with allele frequencies.
