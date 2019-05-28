#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

#Script to extract only the tRNAs that are in the high confidence set from GtRNAdb (as of March 2019)

#####Reads in input files
##############################################################################

my $highconf = $ARGV[0];
my $variants = $ARGV[1];
my $alltRNAs = $ARGV[2];

#####OUTPUT FILES
##############################################################################

system("rm highconf_variants.txt");
open(out1, ">highconf_variants.txt") or die("Cannot open output1 file");
system("rm confidence_tRNAs.txt");
open(out2, ">confidence_tRNAs.txt") or die("Cannot open output2 file");

#####Reads in high confidence tRNA list and stores it
###############################################################################

open(inp0, "$highconf") or die("Cannot open High Confidence tRNA database file");
	my @highconfList;
	
while(<inp0>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq '>'){
		my @line = split(" ", $_);
		
		my @temptRNA = split("_", $line[0]);
		my $tRNA = $temptRNA[2];
		
		push @highconfList, $tRNA;
		
	}
}

close(inp0);

#####Reads in entire tRNA database to annotate high confidence and low confidence
###############################################################################

open(inp2, "$alltRNAs") or die("Cannot open all tRNA database file");
	my %confidencetRNAs;
	my %isoacceptor;
	
while(<inp2>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq '>'){
		my @line = split(" ", $_);
		
		my @temptRNA = split("_", $line[0]);
		my $tRNA = $temptRNA[2];
		
		$confidencetRNAs{$tRNA} = 0;
		
		my @tempIsoacceptor = split("-", $tRNA);
		$isoacceptor{$tRNA} = $tempIsoacceptor[1];
		
		foreach(@highconfList){
			if($tRNA eq $_){
				$confidencetRNAs{$tRNA} = 1;
			}
		}
		
	}
}

close(inp2);

print out2 "tRNA\tConfidence\tIsoacceptor\n";

foreach my $keys (keys %confidencetRNAs){
	if($confidencetRNAs{$keys} == 1){
		print out2 "$keys\tHigh\t$isoacceptor{$keys}\n";
	}
	
	else{
		print out2 "$keys\tLow\t$isoacceptor{$keys}\n";
	}
}


#####Reads in variants, checks if they are in the high confidence list and if they are, outputs them to a file
###############################################################################

open(inp1, "$variants") or die("Cannot open variant tRNA file");

while(<inp1>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq '#'){
		print out1 "$_\n";
	}
	
	if($check ne '#'){
		my $FullLine = $_;
		my @line = split("\t", $_);
		my $tRNA = $line[0];
		
		foreach(@highconfList){
			if($_ eq $tRNA){
				print out1 "$FullLine\n";
			}
		}
	}
}

close(inp1);
