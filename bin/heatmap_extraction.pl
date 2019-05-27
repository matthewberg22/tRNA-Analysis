#!/usr/bin/perl
use strict;
use warnings;

#By: Matthew Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

## This is a script to extract information for each individual how many mutations they have per tRNA (0, 1 or 2)
## Input 1: List of samples (from demographics.txt)
## Input 2: List of all tRNAs (from capturearray_coverage.txt)
## Input 3: allmutants_AF.txt

# Reads input files
my $Samples = $ARGV[0];
my $tRNAList = $ARGV[1];
my $variants = $ARGV[2];

# Output files
system("rm alleleheatmap.txt");
open(out1, ">alleleheatmap.txt") or die("Cannot open output file");
system("rm alleleheatmap_uncommon.txt");
open(out2, ">alleleheatmap_uncommon.txt") or die("Cannot open output file");

# Reads in samples and tRNA lists to initalize matrix (with 0s for everything)
open(inp0, "$Samples") or die("Cannot open sample file");
	my %HM;
	my %HMAF;
	my @samples;
	my @tRNAs;


while(<inp0>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq 'G'){
	
		my @splitLine = split("\t", $_);
		push @samples, $splitLine[0];
	}
}

close(inp0);

open(inp1, "$tRNAList") or die("Cannot open tRNA list");

while(<inp1>){
	chomp;
	
	my $check = substr($_, 0, 5);
	
	if($check eq 'tRNA-'){
	
		my @splitLine = split("\t", $_);
		my @splitLine2 = split(" ", $splitLine[0]);
		push @tRNAs, $splitLine2[0];
	
	}
}

# Initalize all hashes for each individual and all the tRNAs with 0s

foreach my $sample (@samples){
	foreach my $tRNA (@tRNAs){
		
		$HM{$sample}{$tRNA} = 0;
		$HMAF{$sample}{$tRNA} = 0;
		
	}
}



# Saves allele information for each tRNA for each individual
open(inp0, "$variants") or die("Cannot open all variant tRNA file");

while(<inp0>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq 'G'){
		my @splitLine = split("\t", $_);
		my @splitLine2 = split(" ", $splitLine[1]);
		
		my $SampleID = $splitLine[0];
		my $tRNA = substr($splitLine2[0], 1);
		my $counts = $splitLine[4];
		my $AF = $splitLine[5];
		
		$HM{$SampleID}{$tRNA} = $counts;
		
		if($AF <= 0.05){
			print "$AF\n";
			$HMAF{$SampleID}{$tRNA} = $counts;
		}

	
	}
}

foreach(@tRNAs){
	print out1 "\t$_";
}
print out1 "\n";

foreach my $sample (@samples){
	print out1 "$sample";
	foreach my $keys (sort @tRNAs){
		print out1 "\t$HM{$sample}{$keys}";
	}
	print out1 "\n";
}



foreach(@tRNAs){
	print out2 "\t$_";
}
print out2 "\n";

foreach my $sample (@samples){
	print out2 "$sample";
	foreach my $keys (sort @tRNAs){
		print out2 "\t$HMAF{$sample}{$keys}";
	}
	print out2 "\n";
}