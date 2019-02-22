#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#February 8, 2019

#Script to count how many variants for each isoacceptor and isodecoder in the non-redundant variant tRNA set and also counts how many total copies of each gene we see in the sample set
#To run enter: counting_mutant_iso.pl nonredundant_mutants.txt total_counts.txt
#nonredundant_mutants.txt and total_counts.txt come from Allele_Frequencies.pl script

#####Reads in input files
##############################################################################

my $nonredundant = $ARGV[0];
my $totalcounts = $ARGV[1];

#####OUTPUT FILES
##############################################################################
system("rm All_mutants_isoacceptor.txt");
open(out1, ">All_mutants_isoacceptor.txt") or die("Cannot open output1 file");
system("rm All_mutants_isodecoder.txt");
open(out2, ">All_mutants_isodecoder.txt") or die("Cannot open output2 file");
system("rm Uncommon_mutants_isoacceptor.txt");
open(out3, ">Uncommon_mutants_isoacceptor.txt") or die("Cannot open output3 file");
system("rm Uncommon_mutants_isodecoder.txt");
open(out4, ">Uncommon_mutants_isodecoder.txt") or die("Cannot open output4 file");
system("rm TotalCounts_isoacceptor.txt");
open(out5, ">TotalCounts_isoacceptor.txt") or die("Cannot open output5 file");
system("rm TotalCounts_isodecoder.txt");
open(out6, ">TotalCounts_isodecoder.txt") or die("Cannot open output6 file"); 
 
#####Reads in nonredundant tRNA variant list and adds counts to isoacceptor and isodecoder families
###############################################################################

open(inp0, "$nonredundant") or die("Cannot open tRNA_sequence file");
	my %Isoacceptor;
	my %Isodecoder;
	my $AF;
	my %UncommonIsoacceptor;
	my %UncommonIsodecoder;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check ne "#"){
		my @splitLine = split("\t", $_);
		my @tRNA = split("-", $splitLine[0]);
		$AF = $splitLine[2];
		
		my $check2 = substr($tRNA[0], 0, 1);
		
		if($check2 eq "t"){
			$Isoacceptor{$tRNA[1]}++;
			$Isodecoder{$tRNA[2]}++;
			
			if($AF < 0.05){
			$UncommonIsoacceptor{$tRNA[1]}++;
			$UncommonIsodecoder{$tRNA[2]}++;
			}
		
		}
		
		else{
		}
	}
}

close(inp0);

print out1 "Isoacceptor\tAllCount\n";

foreach my $keys (keys %Isoacceptor){
	print out1 "$keys\t$Isoacceptor{$keys}\n";
}

print out2 "Isodecoder\tAllCount\n";

foreach my $keys (keys %Isodecoder){
	print out2 "$keys\t$Isodecoder{$keys}\n";
}

print out3 "Isoacceptor\tUncommonCount\n";

foreach my $keys (keys %UncommonIsoacceptor){
	print out3 "$keys\t$UncommonIsoacceptor{$keys}\n";
}

print out4 "Isodecoder\tUncommonCount\n";

foreach my $keys (keys %UncommonIsodecoder){
	print out4 "$keys\t$UncommonIsodecoder{$keys}\n";
}
		
		
open(inp1, "$totalcounts") or die("Cannot open total counts file");
	my %seenIsoacceptor;
	my %seenIsodecoder;

while(<inp1>){
	chomp;
	my @splitLine = split("\t", $_);
	my @tRNA = split("-", $splitLine[0]);
	
	if(scalar @tRNA > 1){
		my $check = substr($tRNA[0], 0, 1);
		
		if($check eq 't'){
			$seenIsoacceptor{$tRNA[1]} = 0;
			$seenIsodecoder{$tRNA[2]} = 0;
		}
		
		else{
			my $name = $tRNA[0]."-".$tRNA[2];
			$seenIsoacceptor{$name} = 0;
			my $name2 = $tRNA[0]."-".$tRNA[3];
			$seenIsodecoder{$name2} = 0;
		
		}
	}
}

open(inp1, "$totalcounts") or die("Cannot open total counts file");

while(<inp1>){
	chomp;
	my @splitLine = split("\t", $_);
	my @tRNA = split("-", $splitLine[0]);
	
	if(scalar @tRNA > 1){
		my $coverage = $splitLine[1];
		my $check = substr($tRNA[0], 0, 1);
		
		if($check eq 't'){
			$seenIsoacceptor{$tRNA[1]} = $seenIsoacceptor{$tRNA[1]} + $coverage;
			$seenIsodecoder{$tRNA[2]} = $seenIsodecoder{$tRNA[2]} + $coverage;
		}
		
		else{
			my $name = $tRNA[0]."-".$tRNA[2];
			$seenIsoacceptor{$name} = $seenIsoacceptor{$name} + $coverage;
			my $name2 = $tRNA[0]."-".$tRNA[3];
			$seenIsodecoder{$name2} = $seenIsodecoder{$name2} + $coverage;
		
		}
	}
}

print out5 "Isoacceptor\tSeenCounts\n";

foreach my $keys (keys %seenIsoacceptor){
	print out5 "$keys\t$seenIsoacceptor{$keys}\n";
}

print out6 "Isodecoder\tSeenCounts\n";

foreach my $keys (keys %seenIsodecoder){
	print out6 "$keys\t$seenIsodecoder{$keys}\n";
}
		
		
		
		
		
		
		
	
