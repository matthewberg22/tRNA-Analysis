#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#February 12, 2019

#Script to analyze file of tRNAs with more than one allele
#To run, enter: problem_analysis.pl problem_loci.txt
#problem_loci.txt comes from tRNA_trimming.pl script

#####Reads in input files
##############################################################################

my $problems = $ARGV[0];

######OUTPUT FILES
##############################################################################

system("rm morethan2Individual.txt"); 
open(out1, ">morethan2Individual.txt") or die("Cannot open output1 file");
system("rm morethan2tRNA.txt");
open(out2, ">morethan2tRNA.txt") or die("Cannot open output2 file");
system("rm morethan2howmany.txt");
open(out3, ">morethan2howmany.txt") or die("Cannot open output2 file");

#####Reads in counts file and sums up totals
###############################################################################

open(inp0, "$problems") or die("Cannot open allele_counts file");
	my %ID;
	my %tRNAcounts;
	my %totalAllele;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 2);
	
	if($check eq 'GL'){
		$ID{$_}++;
	}
	
	elsif($check eq "tR"){
		my @splitLine = split("\t", $_);
		$tRNAcounts{$splitLine[0]}++;
		$totalAllele{$splitLine[1]}++;
	}
}

close(inp0);

foreach my $keys (keys %ID){
	print out1 "$keys\t$ID{$keys}\n";
}

foreach my $keys (keys %tRNAcounts){
	print out2 "$keys\t$tRNAcounts{$keys}\n";
}

foreach my $keys (keys %totalAllele){
	print out3 "$keys\t$totalAllele{$keys}\n";
}

		
