#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

#Script to extract how many tRNAs were identified for each individual and list of tRNAs that are missing or low coverage
#To run enter extract_analysis.pl BLAST_analysis.txt
#BLAST_analysis.txt comes from BLAST_analysis.pl and there will be one file for each sample in dataset

#####Reads in input files
##############################################################################

my $BLAST_analysis = $ARGV[0]; 

#####Reads in BLAST_analysis file and prints out patient ID, # of tRNAs identified and a list of tRNAs missing
###############################################################################

open(inp0, "$BLAST_analysis") or die("Cannot open BLAST_analysis file");
	my $Identifier;
	my $TotaltRNAs;
	my $TotalReads;
	my $ConfidentReads;
	my $AmbiguousReads;
	my $DuplicatedtRNAs;


while(<inp0>){
	chomp;
	my $check = substr($_, 0, 21);
	
	if($check eq "***** IDENTIFER No.: "){
		my @tempIdentifier = split(":", $_);
		my @tempIdentifier2 = split(" ", $tempIdentifier[1]);
		$Identifier = $tempIdentifier2[0];
	}
	
	if($check eq "Total tRNAs identifie"){
		my @tempTotaltRNAs = split(":", $_);
		my @tempTotaltRNAs2 = split(" ", $tempTotaltRNAs[1]);
		$TotaltRNAs = $tempTotaltRNAs2[0];
	}
	
	if($check eq "Total reads with full"){
		my @tempTotalReads = split(":", $_);
		my @tempTotalReads2 = split(" ", $tempTotalReads[1]);
		$TotalReads = $tempTotalReads2[0];
	}
	
	if($check eq "Total confident reads"){
		my @tempConfidentReads = split(":", $_);
		my @tempConfidentReads2 = split(" ", $tempConfidentReads[1]);
		$ConfidentReads = $tempConfidentReads2[0];
	}
	
	if($check eq "Total ambiguous reads"){
		my @tempAmbiguousReads = split(":", $_);
		my @tempAmbiguousReads2 = split(" ", $tempAmbiguousReads[1]);
		$AmbiguousReads = $tempAmbiguousReads2[0];
	}
	
	if($check eq "Reads with two tRNAs:"){
		my @tempDuplicatedtRNAs = split(":", $_);
		my @tempDuplicatedtRNAs2 = split(" ", $tempDuplicatedtRNAs[1]);
		$DuplicatedtRNAs = $tempDuplicatedtRNAs2[0];
	}
}

#Results are printed to the screen so they can be piped into a single file for all the individuals
print "$Identifier\t$TotaltRNAs\t$TotalReads\t$ConfidentReads\t$AmbiguousReads\t$DuplicatedtRNAs\n";
	
	



