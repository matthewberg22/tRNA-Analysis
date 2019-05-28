#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

#Script to analyze file of nonredundant tRNAs mutants
#To run, enter: nonredundant_tRNA_analysis.pl nonredundant_mutants.txt
#nonredundant_mutants.txt comes from Allele_Frequencies.pl script
#Output is: file of all mutations for each tRNA loci (or group of loci if they could not be confidently identified), file of uncommon mutations for each for each tRNA loci (or group of loci if they could not be confidently identified), file of how many mutations occur in each allele

#####Reads in input files
##############################################################################

my $nonredundant = $ARGV[0];

######OUTPUT FILES
##############################################################################
system("rm Allmutation_tRNAcount.txt");
open(out1, ">Allmutation_tRNAcount.txt") or die("Cannot open output1 file");
system("rm Uncommonmutation_tRNAcount.txt");
open(out2, ">Uncommonmutation_tRNAcount.txt") or die("Cannot open output2 file");
system("rm Mutationsperallele.txt");
open(out3, ">Mutationsperallele.txt") or die("Cannot open output2 file");

#####Reads in counts file and sums up totals
###############################################################################

open(inp0, "$nonredundant") or die("Cannot open file");
	my %tRNAname;
	my %tRNAnamecounter;
	my %AlleleFreq;
	my %MutsPertRNA;
	my %Uncommon;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq 't'){
		my @splitLine = split("\t", $_);
		my $tRNA = $splitLine[0];
		$tRNAname{$splitLine[1]} = $tRNA;
		$tRNAnamecounter{$tRNA}++;
		$AlleleFreq{$splitLine[1]} = $splitLine[2];
	}
	
	else{
	}
}

close(inp0);


#Prints the counts of mutations seen for each tRNA
foreach my $keys (keys %tRNAnamecounter){
	print out1 "$keys\t$tRNAnamecounter{$keys}\n";
}

#Prints the counts of mutations seen for each tRNA with allele frequency less than 0.05
foreach my $keys (keys %tRNAname){
	if($AlleleFreq{$keys} < 0.05){
		my $tRNA = $tRNAname{$keys};
		$Uncommon{$tRNA}++;
	}
}

foreach my $keys (keys %Uncommon){
	print out2 "$keys\t$Uncommon{$keys}\n";
}

#Calculates how many alleles have 1 mutation, 2 mutations, 3 mutations, etc.
foreach my $keys (keys %tRNAname){
	my @tRNA = split(" ", $keys);
	my $numMuts = scalar @tRNA;
	$MutsPertRNA{$numMuts}++;
}

#Prints out the counts table of  # mutations and counts of how many alleles
foreach my $keys (sort keys %MutsPertRNA){
	print out3 "$keys\t$MutsPertRNA{$keys}\n";
}

		
