#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#February 8, 2019

# Script to count coverage for each tRNA targeted for sequence
# To run enter script.pl all_coverage.txt

#####Reads in input files
##############################################################################

my $coverage = $ARGV[0];
my $tRNAsequence = $ARGV[0];

#####OUTPUT FILES
##############################################################################
open(out1, ">capturearray_coverage.txt") or die("Cannot open output1 file");


#####Reads in BLAST_analysis file and prints out patient ID, # of tRNAs identified and a list of tRNAs missing
###############################################################################

open(inp0, "$coverage") or die("Cannot open coverage file");
	my %totalcoverage;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq "G"){
		my @splitLine = split("\t", $_);
		
		my $tRNAname = $splitLine[1];
		my $coverage = $splitLine[3];
		
		if($coverage > 9){
		
		$totalcoverage{$tRNAname} = 0;
		
		}
	}
}

open(inp0, "$coverage") or die("Cannot open coverage file");


while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq "G"){
		my @splitLine = split("\t", $_);
		
		my $tRNAname = $splitLine[1];
		my $coverage = $splitLine[3];
		
		if($coverage > 9){
		
		$totalcoverage{$tRNAname} = $totalcoverage{$tRNAname} + $coverage;
		
		}
	}
}

close(inp0);

print out1 "tRNA\tTotalCoverage\n";



open(inp1, "$tRNAsequence") or die("Cannot open tRNA_sequence file");
	my %tRNAgene;

while(<inp1>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq ">"){
		my @splitLine = split("\t", $_);
		my @tRNA = split("_", $splitLine[0]);
		
		$tRNAgene{$tRNA[2]} = 0;
	
	}
}

close(inp1);

foreach my $keys (keys %totalcoverage){
	my @tRNAgroup = split(" ", $keys);
	
	foreach(@tRNAgroup){
		$tRNAgene{$_} = 1;
	}
}

foreach my $keys (keys %tRNAgene){
	if($tRNAgene{$keys} == 0){
		$totalcoverage{$keys} = 0;
	}
}

foreach my $keys (keys %totalcoverage){
	print out1 "$keys\t$totalcoverage{$keys}\n";
}
