#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

#Script to count coverage for each tRNA targeted for sequencing
#To run enter script.pl Total_Coverage.txt GtRNAdb_20bpflanking.txt
#Total_Coverage.txt comes from BLAST_analysis.pl
#GtRNAdb_20bpflanking.txt is the database file used as a query against the sequencing read BLAST database

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

my $coverage = $ARGV[0];
my $tRNAsequence = $ARGV[1];

#####OUTPUT FILES
##############################################################################
system("rm capturearray_coverage.txt");
open(out1, ">capturearray_coverage.txt") or die("Cannot open output1 file");

#####Reads in Total_Coverage file, initalizes hash for each tRNA name and then reads through file again to sum up total coverage
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

#####Reads in GtRNAdb file and stores names of all tRNA loci designed into capture array
###############################################################################
open(inp1, "$tRNAsequence") or die("Cannot open tRNA_sequence file");
	my %tRNAgene;

while(<inp1>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq ">"){
		my @splitLine = split("\t", $_);
		my @tRNA = split("_", $splitLine[0]);
		my @tRNA2 = split(" ", $tRNA[2]);
		
		$tRNAgene{$tRNA2[0]} = 0;
		
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
	my @multipletRNAs = split(" ", $keys);
	
	my $division = scalar @multipletRNAs;
	foreach(@multipletRNAs){
		my $correctedcoverage = $totalcoverage{$keys}/$division;
		print out1 "$_\t$correctedcoverage\n";
	}
}
