## R Script to make all the graphs for tRNASeq Paper
## May 2019
## Matthew Berg

# Librarys to make the graphs and clean the data
library(ggplot2)
library(gplots)
library(tidyr)
library(cowplot)
library(ggrepel)
library(ggsci)
library(scales)
library(compositions)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)


# Set the directory to the top directory of the Github
setwd()

################################################
#### Plot read merging and BLAST statistics ####
################################################

d <- read.table("usearch-output/stats.txt", sep = "\t", row.names= 1, header = TRUE)
d$Sample <- rownames(d)

# Calculate proportion of reads that passed merging
d$prop <- (d$merged_reads / d$total_reads) * 100
d <- d[d$prop >50,]

# Graphs the total reads per sample
TotalReads <- ggplot(d, aes(Sample, (total_reads*2))) + geom_point(col = "dodgerblue3", size = 2) + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + xlab('') + ylab('Total Reads')  + coord_cartesian(ylim = c(1500000, 5500000))

# Graphs the percent merged
PercentMerged <- ggplot(d, aes(Sample, prop)) + geom_point(col = "dodgerblue3", size = 2) + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + xlab('') + ylab('Percent Merged') + coord_cartesian(ylim = c(50,100))

# Plots total number of tRNAs identified for each sample
data <- read.table("data/summarized_output_Cov10.txt", header=T, sep = "\t")
data$Sample <- factor(data$Identifier, levels = data$Identifier[order(data$Group, data$Identifier)])

dataMerge <- merge(data, d, by = "Sample")

TotaltRNAs <- ggplot(dataMerge, aes(Sample, TotaltRNAs)) + geom_point(size = 2, col = "dodgerblue3") + xlab('Individuals') + ylab('# tRNAs Identified') + theme(axis.text.x=element_blank(), axis.ticks=element_blank())

# Plots total number of reads containing a full length tRNA plus flanking for each sample
GoodReads <- ggplot(dataMerge, aes(Sample,TotalReads)) + geom_point(size = 2, col = "dodgerblue3") + xlab('Individuals') + ylab('# of Reads Containing\n Full Length tRNA') + theme(axis.text.x=element_blank(), axis.ticks=element_blank())

# Export Figure S1 svg
svg(file = "FigureS1.svg", width = 10, height = 8)
plot_grid(TotalReads, PercentMerged, TotaltRNAs, GoodReads, labels = c("A", "B", "C", "D"), align = 'v')
dev.off()

# Summary numbers for the read merging and BLAST analysis mapping
MergeMean <- mean(dataMerge$prop)
MergeSD <- sd(dataMerge$prop)
tRNAMean <- mean(dataMerge$TotaltRNAs)
tRNASD <- sd(dataMerge$TotaltRNAs)
FullReadsMean <- mean(dataMerge$TotalReads)
FullReadsSD <- sd(dataMerge$TotalReads)

################################################################
#### Plotting Distribution of Coverages to determine cutoff ####
################################################################

# Reads in coverage summary file
all.coverage <- read.table("data/blast_analysis_output/TotalCoverage.txt", header = F, sep = "\t")
names(all.coverage) <- c("ID", "tRNA", "Allele", "Coverage")

# Reformats data
formatted <- dcast(all.coverage, tRNA + Allele ~ ID)
rows <- paste(formatted$tRNA, formatted$Allele, sep = "-")
rownames(formatted) <- rows
formatted.final <- formatted[,-2]
formatted.final <- formatted.final[,-1]
formatted.final[is.na(formatted.final)] <- 0
x <- as.vector(clr(formatted.final))

# Plots the frequency of coverage for all unique sequences seen in the 84 individuals sequenced
# Limits the x-axis to 300 to scale plot better
# Red line indicates 10x coverage cut off, blue is 15x and orange is 5x
Allcoverageplot <- ggplot(all.coverage, aes(x = Coverage)) + 
  geom_histogram(binwidth = 1, color = "black", fill = "white") + 
  coord_cartesian(xlim=c(0, 300), ylim=c(0,1000)) + xlab("Coverage of tRNA allele") + 
  ylab("Count") + 
  geom_vline(aes(xintercept = 5), color = 'orange') + 
  geom_vline(aes(xintercept = 10), color = 'red') + 
  geom_vline(aes(xintercept = 15), color = 'blue') +
  theme(legend.position="none")

# Centered log ratio transformation of the data and plot
all.coverage$CLR <- as.vector(clr(all.coverage$Coverage))
clr <- all.coverage[which(all.coverage$Coverage == 5), 5]
clr5 <- clr[1]
clr <- all.coverage[which(all.coverage$Coverage == 10), 5]
clr10 <- clr[1]
clr <- all.coverage[which(all.coverage$Coverage == 15), 5]
clr15 <- clr[1]

CLRcoverageplot <- ggplot(all.coverage, aes(x = CLR)) + 
  geom_histogram(binwidth = 0.5, color = "black", fill = "white") + 
  coord_cartesian(ylim=c(0,30000)) + xlab("clr(Coverage of tRNA allele)") + 
  ylab("Count") + 
  theme(legend.position="none") +
  geom_vline(aes_(xintercept = clr15), color = 'blue') +
  geom_vline(aes_(xintercept = clr10), color = 'red') +
  geom_vline(aes_(xintercept = clr5), color = 'orange')

#Export Figure S2 svg
svg(file = "FigureS2.svg", width = 10, height = 4)
plot_grid(Allcoverageplot, CLRcoverageplot, labels = c("A", "B"))
dev.off()

#################################################
#### Plot Coverage Counts for each tRNA loci ####
#################################################

# Imports data from script output for coverage at each loci and divides each by 84 (for all samples sequenced)
locicoverage <- read.table("data/capturearray_coverage.txt", header = T, sep = "\t")
locicoverage$Coverage.Persample <- locicoverage$TotalCoverage/84

# Plots coverage of each locus as a point
# Note: Groups of loci where reads could not be uniquely mapped were divided by the total number
# of loci grouped together in the PERL script
locicoverage.plot <- ggplot(locicoverage, aes(tRNA, Coverage.Persample)) + geom_point() + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + ylim(0, 500) + geom_hline(aes(yintercept = 10), color = 'red') + xlab("tRNA Loci") + ylab("Coverage per sample")

# Export Figure S3 svg
svg(file = "FigureS3.svg", width = 5, height = 4)
locicoverage.plot
dev.off()

# Summary numbers for the coverage at each locus
MedianCoverage <- median(locicoverage$Coverage.Persample)
MeanCoverage <- mean(locicoverage$Coverage.Persample)
SDCoverage <- sd(locicoverage$Coverage.Persample)

#################################################
#### Plot Distribution of Allele Frequencies ####
#################################################

# Imports list of all variants and their allele frequencies output data from script
AF <- read.table("data/nonredundant_mutants.txt", header = T, sep = "\t")

# Groups variants into sets based on their allele frequency
Uncommon <- AF %>% group_by(AF < 0.05) %>% tally()
Uncommon <- as.numeric(Uncommon[2,2])

Common <- AF %>% group_by(AF < 0.5) %>% tally()
Common <- as.numeric(Common[2,2])
Common <- Common - Uncommon

NotVariants <- AF %>% group_by(AF > 0.5) %>% tally()
NotVariants <- as.numeric(NotVariants[2,2])

# Creates a dataframe of the number of variants in each bin
AF.counts <- data.frame(Uncommon, Common, NotVariants)
AF.counts <- melt(AF.counts)

# Theme needed for the pie chart
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

# Pie chart of the allele frequency distribution for all the variants observed
all <- ggplot(AF.counts, aes(x="", y = value, fill = variable)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
  scale_fill_brewer(name = "Allele Frequencies",palette = "Blues", direction = -1, labels = c("AF < 5%", "5% < AF < 50%", "AF > 50%" )) +
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(x = c(1, 1, 1.2), label = paste0(round(value/cumsum(value)[length(value)]*100), "%"), fontface = "bold"), position = position_stack(vjust = 0.5))

# Imports only the variants in high confidence tRNAs
highconf.AF <- read.table("data/highconf_variants.txt", header = F, sep = "\t")

# Groups variants in the high confidence tRNA set as above
HC.Uncommon <- highconf.AF %>% group_by(V3 < 0.05) %>% tally()
HC.Uncommon <- as.numeric(HC.Uncommon[2,2])

HC.Common <- highconf.AF %>% group_by(V3 < 0.5) %>% tally()
HC.Common <- as.numeric(HC.Common[2,2])
HC.Common <- HC.Common - HC.Uncommon

HC.NotVariants <- highconf.AF %>% group_by(V3 > 0.5) %>% tally()
HC.NotVariants <- as.numeric(HC.NotVariants[2,2])

highconf.AF.counts <- data.frame(HC.Uncommon, HC.Common, HC.NotVariants)
highconf.AF.counts <- melt(highconf.AF.counts)

Highconf <- ggplot(highconf.AF.counts, aes(x="", y = value, fill = variable)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
  scale_fill_brewer(name = "Allele Frequencies",palette = "Blues", direction = -1, labels = c("AF < 5%", "5% < AF < 50%", "AF > 50%" )) +
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(x = c(1, 1, 1.2), label = paste0(round(value/cumsum(value)[length(value)]*100), "%"), fontface = "bold"), position = position_stack(vjust = 0.5))

# Export Figure S4 svg
svg(file = "FigureS4.svg", width = 10, height = 4)
plot_grid(all, Highconf, labels = "AUTO")
dev.off()

##############################
#### Mutations Per Allele ####
##############################

# Imports data from scripts about how many nucleotide changes are in each variant
MutsPerAllele <- read.table("data/Mutationsperallele.txt", header = F, sep = "\t")
names(MutsPerAllele) <- c("NoMutations", "Count")
MutsPerAllele$NoMutations <- factor(MutsPerAllele$NoMutations, levels = c("1", "2", "3", "4"))

MutsPerAllele.plot <- ggplot(MutsPerAllele, aes(x = "", y = Count, fill = NoMutations)) +
  geom_bar(stat = "identity") + 
  coord_polar("y") +
  scale_fill_brewer(name = "Number of\nMutations",palette = "Blues", direction = -1) +
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(x = c(1, 1, 1.2, 1.4), label = paste0(round(Count/cumsum(Count)[length(Count)]*100, digits = 1), "%"), fontface = "bold"), position = position_stack(vjust = 0.5))


##################################
#### tRNA Variants per Person ####
##################################

# Imports the data containing all the variants that each sample had and the demographics file describing each of the samples
data <- read.table("data/allmutants_AF.txt", header = TRUE, sep = "\t")
groups <- read.table("data/demographics.txt", header = TRUE, sep = "\t")

#Assigns variant a count of 1 (to count how many unique deviations from the reference each sample has)
data$counts <- 1

# Sums up how many counts/variants there are per sample and merges it with the demographic information
count.data <- aggregate(counts ~ Sample, data, sum)
count.data <- merge(count.data, groups, by = "Sample")

# Determines how many variants are uncommon (less than 5% allele frequency) and sums how many uncommon variants each sample has, merges it with the dataframe above
uncommon <- data[data$AF <0.05,]
uncommon.aggregate <- aggregate(counts ~ Sample, uncommon, sum)
all.data.counts <- merge(count.data, uncommon.aggregate, by = "Sample")
names(all.data.counts) <- c("Sample", "counts", "Group", "Sex", "uncommon.counts")

# Determines how many variants are less than 25% and sums how many uncommon variants each sample has, merges it with the dataframe above
data.twentyfive <- data[data$AF < 0.25, ]
data.twentyfive.aggregate <- aggregate(counts ~ Sample, data.twentyfive, sum)
all.data.counts <- merge(all.data.counts, data.twentyfive.aggregate, by = "Sample")
names(all.data.counts) <- c("Sample", "counts", "Group", "Sex", "uncommon.counts", "twentyfivepercent.counts")

# Plots the number of variants per person at different AF cut offs
testAF <- seq(from = 0.01, to = 1, by = 0.01)

allelefreqs <- ""
meansmutations <- ""

for(i in testAF){
  test <- data[data$AF < i,]
  test.aggregate <- aggregate(counts ~ Sample, test, sum)
  variantshere <- mean(test.aggregate$counts)
  allelefreqs <- c(allelefreqs, i)
  meansmutations <- c(meansmutations, variantshere)
  i <- i + 0.01
}

allelefreqs <- as.numeric(allelefreqs)
allelefreqs <- allelefreqs[-1]

meansmutations <- as.numeric(meansmutations)
meansmutations <- meansmutations[-1]

mutationsvsAF <- data.frame(allelefreqs, meansmutations)

# Plots the number of variants per person at different allele frequency cut offs
# Red line represents an allele frequency of 25% and the blue line is uncommon variants per person (less than 5%)
mutperAF.plot <- ggplot(mutationsvsAF, aes(x = allelefreqs, y = meansmutations)) + geom_line() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65)) + xlab("Allele Frequency in Sample Set") + ylab("Average Variants Per Individual") + geom_vline(aes(xintercept = 0.05, color = "red")) + theme(legend.position = "none") + geom_vline(aes(xintercept = 0.25, color = "blue"))

# Plots the allele frequency for all the variants in order of lowest allele frequency to highest
# Red line represents an allele frequency of 25% and the blue line is uncommon variants per person (less than 5%)
AFpervariant.plot <- ggplot(data, aes(x = reorder(Mutation, AF), y = AF)) + geom_point() + geom_hline(aes(yintercept = 0.25, color = "blue")) + geom_hline(aes(yintercept = 0.05, color = "red")) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + xlab("Variant") + ylab("Allele Frequency")

# Export Figure S5 svg
svg("FigureS5.svg", width = 10, height = 5)
plot_grid(mutperAF.plot, AFpervariant.plot, labels = "AUTO")
dev.off()

## FIGURE 3 - Part A, B and C ##

# Dotplots of how many variants each individual has at different cut offs, grouped by sex of sample
# Cut off of 25%
Twentyfive.Mutations.sex <- ggplot(data = all.data.counts, aes(y = twentyfivepercent.counts, x = Sex)) + geom_dotplot(aes(fill = Sex), dotsize = 2, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of variant tRNA alleles") + scale_fill_npg() + scale_y_continuous(expand = c(0,0), limits = c(0,45)) + expand_limits(y = 0) + guides(fill = FALSE)

# No cut off (all variants)
AllMutations.Sex <- ggplot(data = all.data.counts, aes(y = counts, x = Sex)) + geom_dotplot(aes(fill = Sex), dotsize = 2, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of total variant tRNA alleles") + scale_fill_npg() + guides(fill = FALSE)

# Uncommon (< 5%) variants)
UncommonMutations.sex <- ggplot(data = all.data.counts, aes(y = uncommon.counts, x = Sex)) + geom_dotplot(aes(fill = Sex), dotsize = 1, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of uncommon (<5%)\nvariant tRNA alleles") + scale_fill_npg() + scale_y_continuous(expand = c(0,0)) + expand_limits(y = 0) + guides(fill = FALSE)

# Export Figure 3 A, B and C svg
svg(file = "Figure3_Part1.svg", width = 12, height = 4)
plot_grid(AllMutations.Sex, Twentyfive.Mutations.sex, UncommonMutations.sex, labels = "AUTO")
dev.off()

## Figure S6 A, B and C ##

#Same as above, but samples grouped by HTG or Control
UncommonMutations.Group <- ggplot(data = all.data.counts, aes(y = uncommon.counts, x = Group)) + geom_dotplot(aes(fill = Group), dotsize = 1, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of uncommon (<5%) variant tRNA alleles") + scale_fill_npg()  + scale_y_continuous(expand = c(0,0)) + expand_limits(y = 0) + guides(fill = FALSE)

Twentyfive.Mutations.Group <- ggplot(data = all.data.counts, aes(y = twentyfivepercent.counts, x = Group)) + geom_dotplot(aes(fill = Group), dotsize = 1, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of common (<25%) variant tRNA alleles") + scale_fill_npg()  + scale_y_continuous(expand = c(0,0)) + expand_limits(y = 0) + guides(fill = FALSE)

AllMutations.Group <- ggplot(data = all.data.counts, aes(y = counts, x = Group)) + geom_dotplot(aes(fill = Group), dotsize = 1, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of total variant tRNA alleles") + scale_fill_npg()  + scale_y_continuous(expand = c(0,0)) + expand_limits(y = 0) + guides(fill = FALSE)

# Export Figure S6 A, B and C svg
svg(file = "FigureS6.svg", width = 12, height = 4)
plot_grid(AllMutations.Group, Twentyfive.Mutations.Group, UncommonMutations.Group, labels = "AUTO")
dev.off()

# Summary numbers of variants per person
PerPersonMean <- mean(all.data.counts$counts)
PerPersonSD <- sd(all.data.counts$counts)
CommonPerPersonMean <- mean(all.data.counts$twentyfivepercent.counts)
CommonPerPersonSD <- sd(all.data.counts$twentyfivepercent.counts)
UncommonPerPersonMean <- mean(all.data.counts$uncommon.counts)
UncommonPerPersonSD <- sd(all.data.counts$uncommon.counts)

##############################################
#### Heatmap of tRNA mutations per person ####
##############################################

# NOTE: tRNA locus or groups of loci where counts of variants for all samples were zeros were manually removed from the heatmap matrix generated in the perl script "heatmap_extraction.pl" to create alleleheatmap_nozeros.txt
# In total, 338 tRNA locus or groups of loci were removed with zero counts across all samples

heatmap <- read.table("alleleheatmap_nozeros.txt", header = T, sep = "\t")
demographics <- read.table("data/demographics.txt", header = T, sep = "\t")
confidence <- read.table("data/confidence_tRNAs.txt", header = T, sep = "\t")

# Extracts all the tRNA sequences that are annotated as high confidence
high <- confidence[confidence$Confidence == "High", 1]
high <- as.vector(high)
high <- gsub("-", "\\.", high)

# Extracts all the samples that are male
male <- demographics[demographics$Sex == "M", 1]
male <- as.vector(male)

# Extracts all the samples that are female
female <- demographics[demographics$Sex == "F", 1]
female <- as.vector(female)

col.order <- c(male, female)

# Converts the heatmap data to a matrix and orders males on the left of the heatmap and females on the right
rownames(heatmap) <- heatmap[,1]
heatmap[,1] <- NULL
matrix.heatmap <- as.matrix(heatmap)
matrix.heatmap <- t(matrix.heatmap)

matrix.heatmap <- matrix.heatmap[,col.order]

# Sets a color panel for the heatmap, where grey means no mutation, blue means one allele (out of two for a diploid organism) is variant and red means both allele at that locus are variant
mycol <- colorpanel(3, "grey95", "#3C5488FF", "#DC0000FF")

#Orders the columns so males are on the left and females are on the right of the heatmap
col.order <- c(male, female)
matrix.heatmap <- matrix.heatmap[,col.order]

# Labeling of the samples: Males are blue, females are red
cols <- rep("#E64B35FF", ncol(matrix.heatmap))
cols[colnames(matrix.heatmap) %in% male] <- "#4DBBD5FF"
cols[colnames(matrix.heatmap) %in% female] <- "#f73b11ff"

# Labeling of the tRNAs on the heat map: High confidence tRNAs are orange, low confidence tRNAs are green
rows <- rep("#429e88f0", nrow(matrix.heatmap))
rows[rownames(matrix.heatmap) %in% high] <- "#eb5a2cff"

# Export Figure 3D pdf
pdf(file = "Figure3D.pdf", width = 9, height = 6)
heatmap.2(matrix.heatmap, trace = "none", dendrogram = "row", Colv = FALSE, col = mycol,  density.info = "none", ylab = "tRNAs", labRow = FALSE, labCol = FALSE, key.xlab = "No. of\nMutations", key.title = "", key.ylab = "", margins = c(2,1.5), keysize = 0.25, lhei = c(2,5), lwid = c(2,5), sepwidth = c(0,0), ColSideColors = cols, RowSideColors = rows)
dev.off()

# Figure S6D seperating heatmap by HTG and Controls
Control <- demographics[demographics$Group == "Control", 1]
Control <- as.vector(Control)

HTG <- demographics[demographics$Group == "HTG", 1]
HTG <- as.vector(HTG)

col.order.group <- c(Control, HTG)

matrix.heatmap <- matrix.heatmap[,col.order.group]

# Label heatmap by sample: Controls are blue, HTG are red
cols <- rep("#E64B35FF", ncol(matrix.heatmap))
cols[colnames(matrix.heatmap) %in% Control] <- "#4DBBD5FF"
cols[colnames(matrix.heatmap) %in% HTG] <- "#f73b11ff"

# Export FigureS6D pdf
pdf(file = "FigureS6D.pdf", width = 9, height = 6)
heatmap.2(matrix.heatmap, trace = "none", dendrogram = "row", Colv = FALSE, col = mycol,  density.info = "none", ylab = "tRNAs", labRow = FALSE, labCol = FALSE, key.xlab = "No. of\nMutations", key.title = "", key.ylab = "", margins = c(2,1.5), keysize = 0.25, lhei = c(2,5), lwid = c(2,5), sepwidth = c(0,0), ColSideColors = cols, RowSideColors = rows)
dev.off()

## Figure 4 - Check invariant vs variant loci ##

# Import whole heatmap data (including zeros)
heatmapall <- read.table("alleleheatmap.txt", header = T, sep = "\t")

rownames(heatmapall) <- heatmapall$X
heatmapall$X <- NULL
heatmapall.matrix <- as.matrix(heatmapall)

# Number of unmutated genes per isoacceptor
DisttRNA <- as.data.frame(colSums(heatmapall))
colnames(DisttRNA) <- "Mutations"

NoMuttRNA <- DisttRNA %>% group_by(Mutations == 0) %>% tally()
NoMuttRNA <- as.numeric(NoMuttRNA[2,2])

MuttRNA <- DisttRNA %>% group_by(Mutations > 0) %>% tally()
MuttRNA <- as.numeric(MuttRNA[2,2])

MuttRNACounts <- data.frame(NoMuttRNA, MuttRNA)
MuttRNACounts <- melt(MuttRNACounts)

isoacceptor <- data.frame(confidence$tRNA, confidence$Isoacceptor, confidence$Confidence)
isoacceptor$confidence.tRNA <- gsub("-", "\\.", isoacceptor$confidence.tRNA)
colnames(isoacceptor) <- c("tRNA", "Isoacceptor", "Confidence")

DisttRNA$tRNA <- rownames(DisttRNA)
ByIsoacceptor <- merge(DisttRNA, isoacceptor, by = "tRNA")

NoMutsIsoacceptor <- ByIsoacceptor[ByIsoacceptor$Mutations == 0,]
NoMutsIsoacceptor.Dist <- count(NoMutsIsoacceptor, Isoacceptor)
NoMutsConfidence.Dist <- count(NoMutsIsoacceptor, Confidence)

MutsIsoacceptor <- ByIsoacceptor[ByIsoacceptor$Mutations > 0,]
MutsIsoacceptor.Dist <- count(MutsIsoacceptor, Isoacceptor)
MutsConfidence.Dist <- count(MutsIsoacceptor, Confidence)

IsoacceptorMut <- merge(NoMutsIsoacceptor.Dist, MutsIsoacceptor.Dist, by = "Isoacceptor", all.x = T, all.y = T)
colnames(IsoacceptorMut) <- c("Isoacceptor", "NoMuts", "Muts")
IsoacceptorMut[is.na(IsoacceptorMut)] <- 0


# Prop.test to see if the ratio of no mutants to mutants is the same for all the isoacceptors
stats.test <- IsoacceptorMut
stats.test$Total <- stats.test$NoMuts + stats.test$Muts
prop.test(stats.test$Muts, stats.test$Total)

# Plot graph of variant vs invariant as stacked bar chart for each isoacceptor
IsoacceptorMut <- melt(IsoacceptorMut)

ConfidenceMut <- merge(NoMutsConfidence.Dist, MutsConfidence.Dist, by = "Confidence")
colnames(ConfidenceMut) <- c("Confidence", "NoMuts", "Muts")
ConfidenceMut <- melt(ConfidenceMut)

MutationPlot <- ggplot(IsoacceptorMut, aes(x = reorder(Isoacceptor, -value), y = value, fill = variable)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + xlab("tRNA Isoacceptor") + ylab("# of genes") + scale_fill_manual(values = c("black", "red"), name = "", labels = c("No Mutations", "Mutations"))

ConfidencePlot <- ggplot(ConfidenceMut, aes(x = reorder(Confidence, -value), y = value, fill = variable)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + xlab("tRNA Isoacceptor") + ylab("# of genes") + scale_fill_manual(values = c("black", "red"), name = "", labels = c("No Mutations", "Mutations"))

gcircle <- grid::circleGrob()
circ <- ggdraw(gcircle)

# Export Figure4 A, B svg
svg(file = "Fig4AB.svg", width = 9, height = 5)
bottomrow <- plot_grid(ConfidencePlot, circ, labels = c("B", "C"), rel_widths = c(1, 0.5))
plot_grid(MutationPlot, bottomrow, labels = c("A", ""), ncol = 1)
dev.off()

##ACCUMULATION PLOT
library(vegan)

sp <- specaccum(heatmapall.matrix, method = "random", permutations = 100)

# Export Figure 4C pdf
pdf("fig4C.pdf", width = 4, height = 4)
plot(sp, log = "xy", ci.type = "poly", col = "red", lwd = 2, ci.lty = 0, ci.col = "lightgrey", xlab = "Number of Individuals", ylab = "tRNAs with mutations", ylim = c(40, 500))
dev.off()

#####################################
#### Plots tRNA Position HeatMap ####
#####################################

highconf <- read.table("data/point_one_map.txt", header = T, sep = "\t")
highconf$Position <- factor(highconf$Position, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "20A", "20B", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "41A", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73"))
highconf.melt <- melt(highconf)

highconf.linegraph <- ggplot(highconf, aes(x = Position)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Count of mutations", color = '')  + geom_line(aes(y = AllCounts, color = "a"), group = 1, size = 1) + geom_line(aes(y = UncommonCounts, color = "b"), linetype = "dashed", group = 1, size = 1) + scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), expand = c(0,0)) + theme(axis.text=element_text(size=8)) + scale_color_manual(name = "Colors", values = c("a" = "black", "b" = "red"))

highconf.linegraph <- ggplot(highconf.melt, aes(Position, variable)) + geom_tile(aes(fill = value), color = "white") + scale_fill_gradientn(colors = c("blue", "white", "white", "white", "red"), values = rescale(c(0, 2.5, 3.5, 4.5, 7))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Export Figure 5 parts svg file
svg(file = "Figure5parts.svg", width = 12, height = 4)
plot_grid(highconf.heatmap, highconf.linegraph, labels = "AUTO", nrow = 2)
dev.off()

##################################
#### Plots scan score changes ####
##################################

Infernal <- read.table("data/infernal_scores.txt", header = T, sep = "\t")
Eufind <- read.table("data/eufind_scores.txt", header = T, sep = "\t")

Infernalscore <- ggplot(Infernal, aes(reorder(Mutation, Ref_Infernal))) + geom_point(aes(y = Ref_Infernal, color = "Reference tRNA"), color = "black") + geom_point(aes(y = Var_Infernal, color = "Variant tRNA"), shape = 21, color = "red") + xlab("tRNA Variant") + ylab("Infernal Score")+ theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + labs(color = "TEST") 

Eufindscore <- ggplot(Eufind, aes(reorder(Mutation, Ref_Eufinder))) + geom_point(aes(y = Ref_Eufinder, color = "Reference tRNA"), color = "black") + geom_point(aes(y = Var_Eufinder, color = "Variant tRNA"), color = "red", shape = 21)   + xlab("tRNA Variant") + ylab("Eufind Score") + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + labs(color = "", alpha = "") + guides(alpha = FALSE)

# Export Figure 6 svg file
svg("Figure6.svg", width = 14, height = 4)
plot_grid(Eufindscore, Infernalscore, labels = c("A", "B"), ncol = 1)
dev.off()
