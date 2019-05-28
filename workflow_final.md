# Workflow for processing tRNA reads
Daniel Giguere and Matthew Berg
Last updated Februrary 8 2019.

**Summary:**

Reads were checked for quality with fastqc. The first 10 and last 2 bases were originally trimmed due to unexpected variations in per base sequence content. This interfered with read overlapping so raw reads were used.

```
# assign variables
FASTQC='/Volumes/data/bin/FastQC/fastqc'
READS='/Volumes/data/trna/reads/'

# output in fastqc_output
$FASTQC $READS -o .
```

# Merge overlapped read pairs with USEARCH

```
usearch='/Volumes/data/bin/usearch'
READS='/Volumes/data/trna/reads'

#test
#$usearch -report usearch_output/GL13.report -fastq_mergepairs $READS/GL13_S1_L001_R1_001.fastq -reverse $READS/GL13_S1_L001_R2_001.fastq -fastqout usearch_output/GL13.fastq -fastaout usearch_output/GL13.fasta


# merge all paired reads with a loop
for f in $READS/*.gz; do

B=`basename $f`
NAME=`echo $B | cut -d "." -f1`

zcat -v $f > $READS/$NAME.fastq

done

for f in $READS/*R1*.fastq; do

#get basename
B=`basename $f`
NAME=`echo $B | cut -d "." -f1`
SAMPLE=`echo $B | cut -d "_" -f1`

R='R2'
# replace R1 with R2 to specify reverse read file
REV=`echo ${B/R1/$R}`
REVNAME=`echo ${NAME/R1/$R}`

$usearch -report usearch_output/$SAMPLE.report -fastq_mergepairs $READS/$NAME.fastq -reverse $READS/$REVNAME.fastq -fastqout usearch_output/$SAMPLE.fastq -fastaout usearch_output/$SAMPLE.fasta

echo $SAMPLE.done

done
```

# Get % of reads merging for each sample

```
#create file
echo sample$'\t'merged_reads$'\t'total_reads > usearch_output/stats.txt

#extract info from each output file
for f in usearch_output/*.report; do

#get sample name
B=`basename $f`
NAME=`echo $B | cut -d "." -f1`

#get merged and total reads number, add to stats file
echo $NAME$'\t'"`sed -n 6p $f | awk '/[0-9]/ {print $1, "\t", $3}'`" >> usearch_output/stats.txt

done

#plot stats in R

R

d <- read.table("usearch_output/stats.txt", sep = "\t", header=TRUE, stringsAsFactors=FALSE)
d$percent <- (d$merged_reads/d$total_reads * 100)

pdf("percent_reads_mapped.pdf")
plot(1:96, d$percent, ylim = c(0,100), main = "Percentage of paired reads merged with USEARCH", xlab = "Sample number", ylab = "Percentage")
dev.off()

#get sample names to be thrown out, write to table
thrown <- d$sample[which(d$percent < 5)]
write.table(as.matrix(thrown), "usearch_output/samples_not_merged.txt", sep = "\n", quote = FALSE, col.names=FALSE, row.names=FALSE)

```

#trouble shoot reads that don't merge


```
#eg GL12471

# get sub sample
$usearch -fastx_subsample GL12471_S6_L001_R1_001.fastq -reverse GL12471_S6_L001_R2_001.fastq -fastqout GL12471_subf.fq -output2 GL12471_subr.fq -sample_size 100

#troubleshooging report
$usearch -fastq_mergepairs GL12471_subf.fq -reverse GL12471_subr.fq -tabbedout out.txt -report report.txt -alnout aln.txt




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

d <- read.table("data/passed_merging.txt", sep = "\t", row.names= 1, header = TRUE)
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
Groups <- read.table("data/Samplegroup.txt"), sep = "\t", header = TRUE)

#Shortens the identifiers to just the GL number
d.t$Sample <- rownames(d.t)
Split <- d.t %>% separate(Sample, c("Sample", NA, NA, NA, NA), "_")
d.t$Sample <- Split$Sample

#Adds which group corresponds to which GL number
d.t$group <- ordered(d.t$Sample, levels = Groups$Sample, labels = Groups$Group)
d.t <- d.t[order(d.t$group),]

#Graphs the percent merged
PercentMerged <- ggplot(d.t, aes(order(d.t$group), prop, col = group, label = Sample)) + geom_point() + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + xlab('') + ylab('Percent Merged') + geom_text_repel(aes(label=ifelse(prop<80,as.character(Sample),'')), force =50) + coord_cartesian(ylim = c(30,100))
svg("percentage_reads_merged_all_samples.svg")
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
for sample in ./*.fasta; do

# get the sample name (i.e., GL141 from GL141_001_R1...)
NAME=`basename $sample | cut -d "." -f1`

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
#setwd in top of analysis directory
all.coverage <- read.table("blast_analysis_output/TotalCoverage.txt", sep = "\t")

colnames(all.coverage) <- c("Patient_id", "tRNA", "allele", "Coverage")
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

# this line has to be run separately.
p1 <- ggplot(all.coverage, aes(x = CLR)) +
  geom_histogram(binwidth = 0.5) +
  coord_cartesian(ylim=c(0,10000)) + xlab("clr(Coverage of tRNA allele)") +
  ylab("Frequency") +
  theme(legend.position="none") +
  geom_vline(aes_(xintercept = clr15), color = 'blue') +
  geom_vline(aes_(xintercept = clr10), color = 'red') +
  geom_vline(aes_(xintercept = clr5), color = 'orange')

# this line has to be run separately.
CoveragePlots <- plot_grid(p0, p1, labels = c("A", "B"))

pdf("CoveragePlots.pdf")
CoveragePlots
dev.off()

```
After looking at the cutoff graphs, cutoff was selected as 10x coverage. Ran perl script again to get list of tRNAs with >10x coverage.

```
nohup ./bin/tRNA_sequences_cutoff10.sh &
```

*NEW*

Ran coverage analysis script to determine coverage of each tRNA (or family of tRNAs) across the capture array.

```
perl capturearray_coverage.pl all_coverage.txt
```

Ran tRNA trimming perl script on all samples to get tRNA mutant sequences without 20 bp 5' flanking, concatenates all individual samples into one file

```

# change tRNA_trimming.pl to the corerct name
./bin/tRNA_trimming_run.sh
# mv output files to data folder
# add sanity check that output exists and looks correct
```

Ran tRNA allele frequency counting script to get count of how many alleles we see for each tRNA total across all individuals and how many times we see each variant tRNA in our sample. Also outputs non-redundant list of mutant tRNA sequences with allele frequencies

```
perl ./bin/Allele_Frequencies.pl total_alleles.txt tRNAmutants.txt
# mv output files to data folder
# add sanity check that output exists and looks correct
# lots of things call from the output, "master file"
```

### Analysis of tRNA Variants

Ran perl script to extract information on how many overall mutations per tRNA, how many uncommon ( AF < 5%) mutations per tRNA and how many mutations per allele.

```
perl ./bin/nonredundant_tRNA_analysis.pl nonredundant_mutants.txt
# mv output files to data folder
# add sanity check that output exists and looks correct
```

### tRNAscan-se

The paths for searching for infernal need to be changed depending on where the install was. Open tRNAscan-SE.conf and change the infernal_dir setting to the location of the cmsearch executable. (https://bioinformatics.stackexchange.com/questions/4637/trnascan-se-error-fatal-unable-to-find-usr-local-bin-cmsearch-executable)

```
infernal_dir: /Volumes/data/bin/infernal/bin
/Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o# -f# -X 1 -L UCSC_tRNA.fasta


/Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o# -f# -X 1 -r# -L new_fasta_nonredundant_mutants.fasta

```


# blast mutants against the reference
makeblast db -type nucl
DB='./all_blast_dbs/GtRNAdb'
nohup blastn -query fasta_nonredundant_mutants-2.fasta -db $DB  -perc_identity 100 -outfmt "7 std btop sseq qlen"  -max_target_seqs 10000 -num_threads 20 > fasta_nonredundant_mutants.blast &

# blast reference against the samples
nohup ./dan_blast.sh &


# run tRNA scan on mutants
nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o# -f# -X 1 -r# -E fasta_nonredundant_mutants.fasta > mutants.eu &

# run tRNA scan on reference using Eufind

nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o# -f# -X 1 -r# -e GtRNAdb.txt &

nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o fasta_nonredund2.eufind -f fasta_nonredund2.fp -X 1 -r# -e fasta_nonredundant_mutants-2.fasta &

# run tRNA scan on reference using Infernal

nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o Gt.infernal -f Gt_struct.infernal -X 1 -r Gt_firstpass.infernal -I GtRNAdb_trimforscore.fasta &

nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o mutants.infernal -f fasta_nonredund2.infernal -X 1 -r fasta_nonredund2_fp.infernal -I -E fasta_nonredundant_mutants-2.fasta &

# run tRNA scan on double mutants using Infernal

nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o dblmutants.infernal -f dblmutant_struct.infernal -X 1 -r dblmutant_firstpass.infernal -I highconf_variants_twomutants.txt &

### Calculate number of reads

for file in ./reads/*.fastq; do

    # print reads, count lines, divide by 4 to get number of reads
    READS=`expr $(cat $file | wc -l) / 4`

    # get sample name
    SAMPLE=`basename $file | cut -d "_" -f1`

    # print to file (tab separated
    echo $SAMPLE$'\t'$READS >> read_counts_per_sample.txt

done


### plotting number of mutations information
```
d <- read.table("alleleheatmap_uncommon.txt", row.names = 1, header = TRUE, sep = "\t", stringsAsFactors=FALSE)

d.sums <- colSums(d)
d.rsums <- rowSums(d)

hist(d.sums, breaks = 30, xlab = "Number of mutations", ylab = "Number of tRNAs", main = "Distribution of mutations among tRNAs")
hist(d.rsums, breaks = 10, xlab = "Number of mutations", ylab = "Number of samples", main = "Distribution of number of mutations among samples")r

````
