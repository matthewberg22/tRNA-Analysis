# Workflow for processing tRNA reads
Daniel Giguere and Matthew Berg

Last updated May 28 2019.

**Summary:**

Reads were checked for quality with fastqc. The first 10 and last 2 bases were originally trimmed due to unexpected variations in per base sequence content. This interfered with read overlapping so raw reads were used.

Download this repository from Github, and set your working directory as the top level directory (i.e., ~/tRNA-analysis). When running these scripts, ensure that your working directory is **always** this top directory for the project.

**Software requirements**

- (Blast suite from NIH)[https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download]. Version 2.2.31 was used for our analysis.
- (USEARCH)[https://www.drive5.com/usearch/]. v11.0.667_i86linux32 was used for our analysis.
- (tRNAscan-SE)[http://lowelab.ucsc.edu/tRNAscan-SE/] Version 2.0 was used for our analysis
-

# Merge overlapped read pairs with USEARCH

```
################################################################################
# Insert the path to your usearch install and downloaded reads.
# usearch='/path/to/usearch'
# READS='/path/to/compressed/reads'
################################################################################

# make directory to store the uncompressed reads. usearch requires uncompressed reads.
mkdir -p reads

# merge all paired reads with a loop
for f in $READS/*.gz; do

B=`basename $f`
NAME=`echo $B | cut -d "." -f1`

# this outputs them in the newly created reads directory
zcat -v $f > $READS/reads/$NAME.fastq

done

# make output directory for usearch files
mkdir -p usearch_output

for f in reads/*R1*.fastq; do

#get basename
B=`basename $f`
NAME=`echo $B | cut -d "." -f1`
SAMPLE=`echo $B | cut -d "_" -f1`

R='R2'
# replace R1 with R2 to specify reverse read file
REV=`echo ${B/R1/$R}`
REVNAME=`echo ${NAME/R1/$R}`

$usearch -report usearch_output/$SAMPLE.report -fastq_mergepairs reads/$NAME.fastq -reverse reads/$REVNAME.fastq -fastqout usearch_output/$SAMPLE.fastq -fastaout usearch_output/$SAMPLE.fasta

echo $SAMPLE.done

done
```

# Get % of reads merging for each sample

```
# create file with headers
echo sample$'\t'merged_reads$'\t'total_reads > usearch_output/stats.txt

# extract info from each output file
for f in usearch_output/*.report; do

# get sample name
B=`basename $f`
NAME=`echo $B | cut -d "." -f1`

# get merged and total reads number, add to stats file
# paste sample name, 6th line contains info we want, print only 1st and 3rd column with awk for line that has numbers
echo $NAME$'\t'"`sed -n 6p $f | awk '/[0-9]/ {print $1, "\t", $3}'`" >> usearch_output/stats.txt

done

```

### Blast against db of reads

For each sample, a BLAST database of the merged reads needs to be made.

Before making blastdb, you need to remove the long string that is common to each of the sample. eg. M00388:388:000000000-BGJ2L from M00388:388:000000000-BGJ2L:1:1101:15326:1913:2. BLAST will cut the identifier to a limited number of characters.

```
# make directory for output
mkdir -p merged_fasta_edited
for sample in usearch_output/*.fasta; do

# get the sample name (i.e., GL141 from GL141_001_R1...)
NAME=`basename $sample | cut -d "." -f1`

# the regex is required since difference samples have different values
sed 's/M00388:[0-9]*:000000000-[a-zA-Z0-9]*//' $sample > merged_fasta_edited/$NAME-merged-edited.fasta

done
```

Make the DBs for each file.

```
# ensure that makeblastdbs.sh is executable
# the rest of the files
# this takes a while to push it to the background.
nohup ./bin/makeblastdbs.sh &
```

Then blast the tRNAs against each database and run the blast anaylsis script.

```
nohup ./bin/blast_unknownCutoff.sh &
```

Next, determine the ideal cutoff by graphing total coverage for every unique tRNA sequence observed from all individuals sequenced. Scripts to make the graphs are found in tRNASeq_Script.R

After looking at the cutoff graphs, cutoff was selected as 10x coverage. Ran perl script again to get list of tRNAs with >10x coverage.

```
nohup ./bin/blast_cutoff10.sh &
```

Ran coverage analysis script to determine coverage of each tRNA (or family of tRNAs) across the capture array. Also, collect summary stats for all 84 samples.

```
perl ./bin/capturearray_coverage.pl ./data/blast_analysis_output/TotalCoverage.txt ./data/GtRNAdb_20flanking.fasta
mv capturearray_coverage.txt ./data/

./bin/extract_summary.sh
```

Run tRNA trimming perl script on all samples to get tRNA mutant sequences without 20 bp 5' flanking or 5 bp 3' flanking, concatenates all individual samples into one file.

```
./bin/tRNA_trimming_run.sh
```

Run tRNA allele frequency counting script to get count of how many alleles we see for each tRNA total across all individuals and how many times we see each variant tRNA in our sample. Also outputs non-redundant list of mutant tRNA sequences with allele frequencies.

```
perl ./bin/Allele_Frequencies.pl ./data/total_alleles.txt ./data/tRNAmutants.txt
mv allmutants_AF.txt ./data/
mv nonredundant_mutants.txt ./data/
mv total_counts.txt ./data/
mv fasta_nonredundant_mutants.fasta ./data/
```

### Analysis of tRNA Variants

Run perl script to extract information on how many overall mutations per tRNA, how many uncommon ( AF < 5%) mutations per tRNA and how many mutations per allele.

```
perl ./bin/nonredundant_tRNA_analysis.pl ./data/nonredundant_mutants.txt
mv Allmutation_tRNAcount.txt ./data/
mv Mutationsperallele.txt ./data/
mv Uncommonmutation_tRNAcount.txt ./data/
```

Convert genomic coordinates into tRNA numbering, based upon standard nomenclature.

```
#Creates a database of all genomic loci for tRNA and the coresponding tRNA number
perl ./bin/individual_tRNA_numberingscript.pl ./data/tRNA_alignment_nonumbers.txt ./data/tRNA_numbering_isoacceptor.txt
mv tRNAstructurenumber.txt ./data/

#Numbers the specific mutants that were identified from our sequencing
perl ./bin/tRNA_variant_annotation.pl ./data/GtRNAdb.txt ./data/tRNAstructurenumber.txt ./data/nonredundant_mutants.txt
mv Annotated_nonredundant_mutants.txt ./data/
```

Subset list of variants to create a new list of only variants occuring in high confidence tRNAs (as defined by GtRNAdb). Run second script to count how many mutations are at each position in the tRNA for all variants in the high confidence set with only one mutation.

```
perl ./bin/highconf_tRNA_variants.pl ./data/highconf_GtRNAdb.fasta ./data/Annotated_nonredundant_mutants.txt ./data/GtRNAdb.txt
mv highconf_variants.txt ./data/
mv confidence_tRNAs.txt ./data/

perl ./bin/mutations_tRNAstructure.pl ./data/highconf_variants.txt
mv point_one_map.txt ./data/
mv indel_one_map.txt ./data/
```

Create a matrix for the heatmap of each individuals tRNA profile.

```
perl ./bin/heatmap_extraction.pl ./data/demographics.txt ./data/capturearray_coverage.txt ./data/allmutants_AF.txt
mv alleleheatmap.txt ./data/
mv alleleheatmap_uncommon.txt ./data/
```

### tRNAscan-se

Run tRNAscan (both Eufind and Infernal seperately) to calculate scores for each reference tRNA and the variant tRNA sequence.

The paths for searching for infernal need to be changed depending on where the install was. Open tRNAscan-SE.conf and change the infernal_dir setting to the location of the cmsearch executable. (https://bioinformatics.stackexchange.com/questions/4637/trnascan-se-error-fatal-unable-to-find-usr-local-bin-cmsearch-executable)

```
infernal_dir: /Volumes/data/bin/infernal/bin
/Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o# -f# -X 1 -L GtRNAdb.fasta


/Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o# -f# -X 1 -r# -L new_fasta_nonredundant_mutants.fasta

```

# Run tRNA scan on reference and mutants using Eufind

```
nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o# -f# -X 1 -r# -e data/GtRNAdb.txt &

nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o fasta_nonredund.eufind -f fasta_nonredund.fp -X 1 -r# -e data/fasta_nonredundant_mutants.fasta &
```

# run tRNA scan on reference using Infernal

```
nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o Gt.infernal -f Gt_struct.infernal -X 1 -r Gt_firstpass.infernal -I data/GtRNAdb_20flanking.fasta &

nohup /Volumes/data/bin/trnascanse/bin/tRNAscan-SE -o mutants.infernal -f fasta_nonredund.infernal -X 1 -r fasta_nonredund_fp.infernal -I -E data/fasta_nonredundant_mutants.fasta &

```
