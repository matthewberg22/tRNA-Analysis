## R Script to make all the graphs for tRNASeq Paper
## May 2019
## Matthew Berg

#Librarys to make the graphs and clean the data
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


#Set the directory to the top directory of the Github
setwd()

################################################
#### Plot read merging and BLAST statistics ####
################################################

d <- read.table("data/stats.txt", sep = "\t", row.names= 1, header = TRUE)
d$Sample <- rownames(d)


# calculate proportion of reads that passed merging
d$prop <- (d$merged_reads / d$total_reads) * 100
d <- d[d$prop >50,]

#Graphs the total reads per sample
TotalReads <- ggplot(d, aes(Sample, (total_reads*2))) + geom_point(col = "dodgerblue3", size = 2) + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + xlab('') + ylab('Total Reads')  + coord_cartesian(ylim = c(1500000, 5500000))


#Graphs the percent merged
PercentMerged <- ggplot(d, aes(Sample, prop)) + geom_point(col = "dodgerblue3", size = 2) + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + xlab('') + ylab('Percent Merged') + coord_cartesian(ylim = c(50,100))

#Plots total number of tRNAs identified for each sample
data <- read.table("data/summarized_output_Cov10.txt", header=T, sep = "\t")
data$Sample <- factor(data$Identifier, levels = data$Identifier[order(data$Group, data$Identifier)])

dataMerge <- merge(data, d, by = "Sample")

TotaltRNAs <- ggplot(dataMerge, aes(Sample, TotaltRNAs)) + geom_point(size = 2, col = "dodgerblue3") + xlab('Individuals') + ylab('# tRNAs Identified') + theme(axis.text.x=element_blank(), axis.ticks=element_blank())

#Plots total number of good reads for each sample
GoodReads <- ggplot(dataMerge, aes(Sample,TotalReads)) + geom_point(size = 2, col = "dodgerblue3") + xlab('Individuals') + ylab('# of Reads Containing\n Full Length tRNA') + theme(axis.text.x=element_blank(), axis.ticks=element_blank())

#Export Figure S1 svg
svg(file = "FigureS1.svg", width = 10, height = 8)
plot_grid(TotalReads, PercentMerged, TotaltRNAs, GoodReads, labels = c("A", "B", "C", "D"), align = 'v')
dev.off()

MergeMean <- mean(dataMerge$prop)
MergeSD <- sd(dataMerge$prop)
tRNAMean <- mean(dataMerge$TotaltRNAs)
tRNASD <- sd(dataMerge$TotaltRNAs)
FullReadsMean <- mean(dataMerge$TotalReads)
FullReadsSD <- sd(dataMerge$TotalReads)

############################################
#### Plotting Distribution of Coverages ####
############################################

#Reads in coverage summary file that is the output from my perl script
all.coverage <- read.table("data/blast_analysis_output/TotalCoverage.txt", header = F, sep = "\t")
names(all.coverage) <- c("ID", "tRNA", "Allele", "Coverage")

formatted <- dcast(all.coverage, tRNA + Allele ~ ID)
rows <- paste(formatted$tRNA, formatted$Allele, sep = "-")
rownames(formatted) <- rows
formatted.final <- formatted[,-2]
formatted.final <- formatted.final[,-1]
formatted.final[is.na(formatted.final)] <- 0
x <- as.vector(clr(formatted.final))


#Plots the frequency of coverage for all unique sequences seen in the 96 individuals sequenced
#Limits the x-axis to 300 to scale plot better
#Red line indicates the mean coverage
#Blue line is my coverage cut off right now
Allcoverageplot <- ggplot(all.coverage, aes(x = Coverage)) + 
  geom_histogram(binwidth = 1, color = "black", fill = "white") + 
  coord_cartesian(xlim=c(0, 300), ylim=c(0,1000)) + xlab("Coverage of tRNA allele") + 
  ylab("Count") + 
  geom_vline(aes(xintercept = 5), color = 'orange') + 
  geom_vline(aes(xintercept = 10), color = 'red') + 
  geom_vline(aes(xintercept = 15), color = 'blue') +
  theme(legend.position="none")

#Centered log ratio transformation of the data
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

#Export svg default
svg(file = "FigureS2.svg", width = 10, height = 4)
plot_grid(Allcoverageplot, CLRcoverageplot, labels = c("A", "B"))
dev.off()

#################################################
#### Plot Coverage Counts for each tRNA loci ####
#################################################

locicoverage <- read.table("data/capturearray_coverage.txt", header = T, sep = "\t")
locicoverage$Coverage.Persample <- locicoverage$TotalCoverage/84

locicoverage.plot <- ggplot(locicoverage, aes(tRNA, Coverage.Persample)) + geom_point() + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + ylim(0, 500) + geom_hline(aes(yintercept = 10), color = 'red') + xlab("tRNA Loci") + ylab("Coverage per sample")

locicoverage.bar <- ggplot(locicoverage, aes(tRNA, Coverage.Persample)) + geom_bar(stat = "identity") + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + ylim(0, 500) + geom_hline(aes(yintercept = 10), color = 'red') + xlab("tRNA Loci") + ylab("Coverage per sample")

#Export svg default
svg(file = "FigureS3.svg", width = 5, height = 4)
locicoverage.plot
dev.off()

MedianCoverage <- median(locicoverage$Coverage.Persample)
MeanCoverage <- mean(locicoverage$Coverage.Persample)
SDCoverage <- sd(locicoverage$Coverage.Persample)


#################################################
#### Plot Distribution of Allele Frequencies ####
#################################################

AF <- read.table("data/nonredundant_mutants.txt", header = T, sep = "\t")

Uncommon <- AF %>% group_by(AF < 0.05) %>% tally()
Uncommon <- as.numeric(Uncommon[2,2])

Common <- AF %>% group_by(AF < 0.5) %>% tally()
Common <- as.numeric(Common[2,2])
Common <- Common - Uncommon

NotVariants <- AF %>% group_by(AF > 0.5) %>% tally()
NotVariants <- as.numeric(NotVariants[2,2])

AF.counts <- data.frame(Uncommon, Common, NotVariants)
AF.counts <- melt(AF.counts)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

all <- ggplot(AF.counts, aes(x="", y = value, fill = variable)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
  scale_fill_brewer(name = "Allele Frequencies",palette = "Blues", direction = -1, labels = c("AF < 5%", "5% < AF < 50%", "AF > 50%" )) +
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(x = c(1, 1, 1.2), label = paste0(round(value/cumsum(value)[length(value)]*100), "%"), fontface = "bold"), position = position_stack(vjust = 0.5))

highconf.AF <- read.table("data/highconf_variants.txt", header = F, sep = "\t")

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

setwd("C:/Users/Matt/Desktop/tRNA Sequencing Manuscript")
svg(file = "FigureS4.svg", width = 10, height = 4)
plot_grid(all, Highconf, labels = "AUTO")
dev.off()

##############################
#### Mutations Per Allele ####
##############################

MutsPerAllele <- read.table("data/Mutationsperallele.txt", header = F, sep = "\t")
names(MutsPerAllele) <- c("NoMutations", "Count")
MutsPerAllele$NoMutations <- factor(MutsPerAllele$NoMutations, levels = c("1", "2", "3", "4"))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

MutsPerAllele.plot <- ggplot(MutsPerAllele, aes(x = "", y = Count, fill = NoMutations)) +
  geom_bar(stat = "identity") + 
  coord_polar("y") +
  scale_fill_brewer(name = "Number of\nMutations",palette = "Blues", direction = -1) +
  blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(x = c(1, 1, 1.2, 1.4), label = paste0(round(Count/cumsum(Count)[length(Count)]*100, digits = 1), "%"), fontface = "bold"), position = position_stack(vjust = 0.5))


svg(file = "FigureS5.svg", width = 5, height = 4)
MutsPerAllele.plot
dev.off()

##################################
#### tRNA Variants per Person ####
##################################

data <- read.table("data/allmutants_AF.txt", header = TRUE, sep = "\t")
groups <- read.table("data/demographics.txt", header = TRUE, sep = "\t")

data$counts <- 1
count.data <- aggregate(counts ~ Sample, data, sum)
count.data <- merge(count.data, groups, by = "Sample")

uncommon <- data[data$AF <0.05,]
uncommon.aggregate <- aggregate(counts ~ Sample, uncommon, sum)
all.data.counts <- merge(count.data, uncommon.aggregate, by = "Sample")
names(all.data.counts) <- c("Sample", "counts", "Group", "Sex", "uncommon.counts")

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

mutperAF.plot <- ggplot(mutationsvsAF, aes(x = allelefreqs, y = meansmutations)) + geom_line() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65)) + xlab("Allele Frequency in Sample Set") + ylab("Average Variants Per Individual") + geom_vline(aes(xintercept = 0.05, color = "red")) + theme(legend.position = "none") + geom_vline(aes(xintercept = 0.25, color = "blue"))

AFpervariant.plot <- ggplot(data, aes(x = reorder(Mutation, AF), y = AF)) + geom_point() + geom_hline(aes(yintercept = 0.25, color = "blue")) + geom_hline(aes(yintercept = 0.05, color = "red")) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank()) + xlab("Variant") + ylab("Allele Frequency")

svg("FigureS6.svg", width = 10, height = 5)
plot_grid(mutperAF.plot, AFpervariant.plot, labels = "AUTO")
dev.off()

#FIGURE 2
Twentyfive.Mutations.sex <- ggplot(data = all.data.counts, aes(y = twentyfivepercent.counts, x = Sex)) + geom_dotplot(aes(fill = Sex), dotsize = 2, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of variant tRNA alleles") + scale_fill_npg() + scale_y_continuous(expand = c(0,0), limits = c(0,45)) + expand_limits(y = 0) + guides(fill = FALSE)

AllMutations.Sex <- ggplot(data = all.data.counts, aes(y = counts, x = Sex)) + geom_dotplot(aes(fill = Sex), dotsize = 2, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of common (<25%) variant tRNA alleles") + scale_fill_npg() + guides(fill = FALSE)

UncommonMutations <- ggplot(data = all.data.counts, aes(y = uncommon.counts, x = Group)) + geom_dotplot(aes(fill = Group), dotsize = 1, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of uncommon (<5%) variant tRNA alleles") + scale_fill_npg()  + scale_y_continuous(expand = c(0,0)) + expand_limits(y = 0) + guides(fill = FALSE)

UncommonMutations.sex <- ggplot(data = all.data.counts, aes(y = uncommon.counts, x = Sex)) + geom_dotplot(aes(fill = Sex), dotsize = 1, binaxis = "y", binwidth = 1, stackdir = "center") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3) + ylab("No. of uncommon (<5%)\nvariant tRNA alleles") + scale_fill_npg() + scale_y_continuous(expand = c(0,0)) + expand_limits(y = 0) + guides(fill = FALSE)

svg(file = "Figure3_Part1.svg", width = 12, height = 4)
plot_grid(Twentyfive.Mutations.sex, UncommonMutations.sex, labels = "AUTO")
dev.off()

PerPersonMean <- mean(all.data.counts$counts)
PerPersonSD <- sd(all.data.counts$counts)
UncommonPerPersonMean <- mean(all.data.counts$uncommon.counts)
UncommonPerPersonSD <- sd(all.data.counts$uncommon.counts)

##############################################
#### Heatmap of tRNA mutations per person ####
##############################################

# NOTE: Zeros were manually removed from the heatmap matrix generated in the perl script "heatmap_extraction.pl"

heatmap <- read.table("alleleheatmap_nozeros.txt", header = T, sep = "\t")
demographics <- read.table("data/demographics.txt", header = T, sep = "\t")
confidence <- read.table("data/confidence_tRNAs.txt", header = T, sep = "\t")

high <- confidence[confidence$Confidence == "High", 1]
high <- as.vector(high)
high <- gsub("-", "\\.", high)

male <- demographics[demographics$Sex == "M", 1]
male <- as.vector(male)

female <- demographics[demographics$Sex == "F", 1]
female <- as.vector(female)

rownames(heatmap) <- heatmap[,1]
heatmap[,1] <- NULL
matrix.heatmap <- as.matrix(heatmap)

mycol <- colorpanel(3, "grey95", "#3C5488FF", "#DC0000FF")

#Males are blue, females are red
cols <- rep("#E64B35FF", nrow(matrix.heatmap))
cols[row.names(matrix.heatmap) %in% male] <- "#4DBBD5FF"
cols[row.names(matrix.heatmap) %in% female] <- "#f73b11ff"

#High confidence tRNAs are orange, low confidence tRNAs are green
rows <- rep("#429e88f0", ncol(matrix.heatmap))
rows[colnames(matrix.heatmap) %in% high] <- "#eb5a2cff"

pdf(file = "heatmap_nozeros.pdf", width = 9, height = 6)
heatmap.2(matrix.heatmap, trace = "none", col = mycol,  density.info = "none", xlab = "tRNAs", ylab = "INDIVIDUALS", labRow = TRUE, labCol = TRUE, key.xlab = "No. of\nMutations", key.title = "", key.ylab = "", margins = c(2,1.5), keysize = 0.25, lhei = c(2,5), lwid = c(2,5), sepwidth = c(0,0), ColSideColors = rows, RowSideColors = cols)
dev.off()

hc.rows <- hclust(dist(matrix.heatmap))
plot(hc.rows)

hc.cols <- hclust(dist(t(matrix.heatmap)))
plot(hc.cols)

##ACCUMULATION PLOT
library(vegan)

heatmapall <- read.table("alleleheatmap.txt", header = T, sep = "\t")

rownames(heatmapall) <- heatmapall$X
heatmapall$X <- NULL
heatmapall.matrix <- as.matrix(heatmapall)

sp <- specaccum(heatmapall.matrix, method = "random", permutations = 100)

pdf("fig4C.pdf", width = 4, height = 4)
plot(sp, log = "xy", ci.type = "poly", col = "red", lwd = 2, ci.lty = 0, ci.col = "lightgrey", xlab = "Number of Individuals", ylab = "tRNAs with mutations", ylim = c(40, 500))
dev.off()



#Number of unmutated genes per isoacceptor
DistPeople <- as.data.frame(rowSums(heatmapall))
colnames(DistPeople) <- "Mutations"
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

colorCount <- length(NoMutsIsoacceptor.Dist$Isoacceptor)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

MutsIsoacceptor <- ByIsoacceptor[ByIsoacceptor$Mutations > 0,]
MutsIsoacceptor.Dist <- count(MutsIsoacceptor, Isoacceptor)
MutsConfidence.Dist <- count(MutsIsoacceptor, Confidence)

IsoacceptorMut <- merge(NoMutsIsoacceptor.Dist, MutsIsoacceptor.Dist, by = "Isoacceptor")
colnames(IsoacceptorMut) <- c("Isoacceptor", "NoMuts", "Muts")

#Prop.test to see if the ratio of no mutants to mutants is the same for all the isoacceptors
stats.test <- IsoacceptorMut
stats.test$Total <- stats.test$NoMuts + stats.test$Muts
prop.test(stats.test$Muts, stats.test$Total)

IsoacceptorMut <- melt(IsoacceptorMut)

ConfidenceMut <- merge(NoMutsConfidence.Dist, MutsConfidence.Dist, by = "Confidence")
colnames(ConfidenceMut) <- c("Confidence", "NoMuts", "Muts")
ConfidenceMut <- melt(ConfidenceMut)

MutationPlot <- ggplot(IsoacceptorMut, aes(x = reorder(Isoacceptor, -value), y = value, fill = variable)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + xlab("tRNA Isoacceptor") + ylab("# of genes") + scale_fill_manual(values = c("black", "red"), name = "", labels = c("No Mutations", "Mutations"))

ConfidencePlot <- ggplot(ConfidenceMut, aes(x = reorder(Confidence, -value), y = value, fill = variable)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + xlab("tRNA Isoacceptor") + ylab("# of genes") + scale_fill_manual(values = c("black", "red"), name = "", labels = c("No Mutations", "Mutations"))

svg(file = "4Aand4B.svg", width = 15, height = 3)
plot_grid(MutationPlot, ConfidencePlot, labels = "AUTO")
dev.off()

#####################################
#### Plots tRNA Position HeatMap ####
#####################################

highconf <- read.table("data/point_one_map.txt", header = T, sep = "\t")
highconf$Position <- factor(highconf$Position, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "20A", "20B", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "41A", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73"))
highconf.melt <- melt(highconf)

highconf.heatmap <- ggplot(highconf.melt, aes(Position, variable)) + geom_tile(aes(fill = value), color = "white") + scale_fill_gsea(breaks = c(0, 2, 4, 6, 8, 10)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

highconf.linegraph <- ggplot(highconf, aes(x = Position)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Count of mutations", color = '')  + geom_line(aes(y = AllCounts, color = "a"), group = 1, size = 1) + geom_line(aes(y = UncommonCounts, color = "b"), linetype = "dashed", group = 1, size = 1) + scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), expand = c(0,0)) + theme(axis.text=element_text(size=8)) + scale_color_manual(name = "Colors", values = c("a" = "black", "b" = "red"))

ggplot(highconf.melt, aes(Position, variable)) + geom_tile(aes(fill = value), color = "white") + scale_fill_gradientn(colors = c("blue", "white", "white", "white", "red"), values = rescale(c(0, 2.5, 3.5, 4.5, 7))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Export svg file
svg(file = "figure2_part1_highconf_crude.svg", width = 12, height = 4)
plot_grid(highconf.heatmap, highconf.linegraph, labels = "AUTO", nrow = 2)
dev.off()

##################################
#### Plots scan score changes ####
##################################

Infernal <- read.table("data/infernal_scores.txt", header = T, sep = "\t")
Eufind <- read.table("data/eufind_scores.txt", header = T, sep = "\t")

Infernalscore <- ggplot(Infernal, aes(reorder(Mutation, Ref_Infernal))) + geom_point(aes(y = Ref_Infernal, color = "Reference tRNA"), color = "black") + geom_point(aes(y = Var_Infernal, color = "Variant tRNA"), shape = 21, color = "red") + xlab("tRNA Variant") + ylab("Infernal Score")+ theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + labs(color = "TEST") 

Eufindscore <- ggplot(Eufind, aes(reorder(Mutation, Ref_Eufinder))) + geom_point(aes(y = Ref_Eufinder, color = "Reference tRNA"), color = "black") + geom_point(aes(y = Var_Eufinder, color = "Variant tRNA"), color = "red", shape = 21)   + xlab("tRNA Variant") + ylab("Eufind Score") + theme(axis.text.x=element_blank(), axis.ticks=element_blank()) + labs(color = "", alpha = "") + guides(alpha = FALSE)

svg("figure2_part2.svg", width = 14, height = 4)
plot_grid(Eufindscore, Infernalscore, labels = c("A", "B"), ncol = 1)
dev.off()
