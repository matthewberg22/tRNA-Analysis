#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#February 6, 2019

# Script to trim count the allele frequency in our sample and the frequency of each mutant tRNA
# To run enter script.pl total_allele_count.txt (from trimming script) mutant_list.txt

#####Reads in input files
##############################################################################

my $counts = $ARGV[0];
my $mutants = $ARGV[1];

#####OUTPUT FILES
##############################################################################

system("rm total_counts.txt"); #for Windows computers use del and for UNIX use rm -f
open(out1, ">total_counts.txt") or die("Cannot open output1 file");
system("rm allmutants_AF.txt"); #for Windows computers use del and for UNIX use rm -f
open(out2, ">allmutants_AF.txt") or die("Cannot open output2 file"); 
system("rm nonredundant_mutants.txt"); #for Windows computers use del and for UNIX use rm -f
open(out3, ">nonredundant_mutants.txt") or die("Cannot open output3 file");  

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

foreach my $keys (sort keys %tRNAtotalcounts){
	print out1 "$keys\t$tRNAtotalcounts{$keys}\n";
}

print "There were $samplecounter individuals in this sample";

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

# foreach my $keys (sort keys %MutantFreq){
	# print "$keys \t $MutantFreq{$keys}\n";
# }

open(inp1, "$mutants") or die("Cannot open allele_counts file");
	my %tRNAID;
	my %AFall;
	my %sequence;
	my $mutation;

while(<inp1>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq "#"){
		print out2 "$_\n";
	}
	
	elsif($check eq ">"){
		my @splitLine = split("\t", $_);
		my $tRNAname = substr($splitLine[0], 1);
		$mutation = $splitLine[1];
		my $AF = $MutantFreq{$mutation}/$tRNAtotalcounts{$tRNAname};
		print out2 "$_\tAF=$AF\n";
		$tRNAID{$mutation} = $tRNAname;
		$AFall{$mutation} = $AF;
	}
	
	elsif($check eq "A"|$check eq "C"|$check eq "G"|$check eq "T"){
		print out2 "$_\n";
		$sequence{$mutation} = $_;
	}
}
close(inp1);

foreach my $keys (sort keys %tRNAID){
	print out3 ">$tRNAID{$keys}\t$keys\tAF=$AFall{$keys}\n";
	print out3 "$sequence{$keys}\n";
}
