#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

#Script to count the allele frequency in our sample and the frequency of each mutant tRNA
#To run enter: Allele_Frequencies.pl total_alleles.txt mutant_list.txt
#Both arguments come from tRNA_trimming.pl

# Copyright (C) 2019 Matthew Berg

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#####Reads in input files
##############################################################################

my $counts = $ARGV[0];
my $mutants = $ARGV[1];

#####OUTPUT FILES
##############################################################################

system("rm total_counts.txt");
open(out1, ">total_counts.txt") or die("Cannot open output1 file");
system("rm allmutants_AF.txt");
open(out2, ">allmutants_AF.txt") or die("Cannot open output2 file"); 
system("rm nonredundant_mutants.txt");
open(out3, ">nonredundant_mutants.txt") or die("Cannot open output3 file");
system("rm fasta_nonredundant_mutants.fasta"); 
open(out4, ">fasta_nonredundant_mutants.fasta") or die("Cannot open output4 file");
  

#####Reads in counts file and sums up totals
###############################################################################

open(inp0, "$counts") or die("Cannot open allele_counts file");
	my %tRNAtotalcounts;
	my $samplecounter = 0;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq '#'){
		$samplecounter++;
	}
	
	if($check ne "#"){
		my @splitLine = split("\t", $_);
		$tRNAtotalcounts{$splitLine[0]} = 0
	}
}

close(inp0);

open(inp0, "$counts") or die("Cannot open allele_counts file");

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check ne "#"){
		my @splitLine = split("\t", $_);
		$tRNAtotalcounts{$splitLine[0]} = $tRNAtotalcounts{$splitLine[0]} + $splitLine[1];
	}
}

close(inp0);

print out1 "tRNA\tTotalAlleles\n";

#Prints the total allele counts for each loci (or group of loci if they cannot be uniquely identified) from all samples
foreach my $keys (sort keys %tRNAtotalcounts){
	print out1 "$keys\t$tRNAtotalcounts{$keys}\n";
}


#####Reads in mutant file and counts up the number of times each mutation occurs
###############################################################################

open(inp1, "$mutants") or die("Cannot open allele_counts file");
	my %MutantFreq;

while(<inp1>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq ">"){
		my @splitLine = split("\t", $_);
		my $tRNAname = $splitLine[1];
		$MutantFreq{$tRNAname} = 0;
	}
}
close(inp1);

open(inp1, "$mutants") or die("Cannot open allele_counts file");

while(<inp1>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq ">"){
		my @splitLine = split("\t", $_);
		my $tRNAname = $splitLine[1];
		my $MutCount = $splitLine[3];
		$MutantFreq{$tRNAname} = $MutantFreq{$tRNAname} + $MutCount;
	}
}
close(inp1);

open(inp1, "$mutants") or die("Cannot open allele_counts file");
	my %tRNAID;
	my %AFall;
	my %sequence;
	my $mutation;

print out2 "Sample\ttRNA\tMutation\tCoverage\tCount\tAF\n";
	
while(<inp1>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq "#"){
		print out2 substr($_, 2);
		print out2 "\t";
	}
	
	elsif($check eq ">"){
		my @splitLine = split("\t", $_);
		my $tRNAname = substr($splitLine[0], 1);
		$mutation = $splitLine[1];
		my $AF = $MutantFreq{$mutation}/$tRNAtotalcounts{$tRNAname};
		print out2 "$_\t$AF\n";
		$tRNAID{$mutation} = $tRNAname;
		$AFall{$mutation} = $AF;
	}
	
	elsif($check eq "A"|$check eq "C"|$check eq "G"|$check eq "T"|$check eq "-"){
		$sequence{$mutation} = $_;
	}
}
close(inp1);

print out3 "Header-tRNA\tMutation\tAF\tSeq\n";

foreach my $keys (sort keys %tRNAID){
	print out3 "$tRNAID{$keys}\t$keys\t$AFall{$keys}\t";
	print out3 "$sequence{$keys}\n";
	
	my @splitkey = split(" ", $keys);
	my @splittRNA = split(" ", $tRNAID{$keys});
	my $singletRNA = $splittRNA[0];
	
	print out4 ">$singletRNA.";
	foreach(@splitkey){
		print out4 "$_.";
	}
	print out4 "\n";
	if(index($sequence{$keys}, "-") != -1){
		$sequence{$keys} =~ tr/-//d;
		}	
		
	print out4 "$sequence{$keys}\n";
}



		
