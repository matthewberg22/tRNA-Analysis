#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#February 8, 2019

#Script to count how many genes for each isoacceptor and isodecoder are present in GtRNAdb file
#To run enter script.pl GtRNAdb.txt

#####Reads in input files
##############################################################################

my $tRNAsequence = $ARGV[0];

#####OUTPUT FILES
##############################################################################
system("rm genes_isoacceptor.txt");
open(out1, ">genes_isoacceptor.txt") or die("Cannot open output1 file");
system("rm genes_isodecoder.txt");
open(out2, ">genes_isodecoder.txt") or die("Cannot open output2 file"); 

#####Reads in GtRNAdb file and adds counts to each isoacceptor or isodecoder family
###############################################################################

open(inp0, "$tRNAsequence") or die("Cannot open tRNA_sequence file");
	my %Isoacceptor;
	my %Isodecoder;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq ">"){
		my @splitLine = split("\t", $_);
		my @tRNA = split("_", $splitLine[0]);
		my @tRNA2 = split("-", $tRNA[2]);
		
		my $check2 = substr($tRNA2[0], 0, 1);
		
		if($check2 eq "t"){
		$Isoacceptor{$tRNA2[1]}++;
		$Isodecoder{$tRNA2[2]}++;
		}
		
		else{
		my $name = $tRNA2[0]."-".$tRNA2[2];
		$Isoacceptor{$name}++;
		my $name2 = $tRNA2[0]."-".$tRNA2[3];
		$Isodecoder{$name2}++;
		}
	}
}

close(inp0);

print out1 "Isoacceptor\tCount\n";

foreach my $keys (keys %Isoacceptor){
	print out1 "$keys\t$Isoacceptor{$keys}\n";
}

print out2 "Isodecoder\tCount\n";

foreach my $keys (keys %Isodecoder){
	print out2 "$keys\t$Isodecoder{$keys}\n";
}

		
		
		
		
		
		
		
		
		
		
	
