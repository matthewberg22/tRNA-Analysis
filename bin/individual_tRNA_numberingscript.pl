#!/usr/bin/perl 
use strict;
use warnings;

#By: Matthew Berg
#Date: May 2019
#Deparment of Biochemistry, University of Western Ontario

#Script to number all tRNA genes in GtRNAdb with canonical tRNA numbering, based off alignments and structure information in GtRNAdb
#Usage is: individual_tRNA_numberingscript.pl tRNA_alignment_nonumbers.txt tRNA_numbering_isoacceptor.txt
#tRNA_numbering_isoacceptor.txt was created manually from the structure/alignment information on GtRNAdb and can be found in the data folder

#Reads in input files
my $structure = $ARGV[0];
my $numbering = $ARGV[1];

#OUTPUT FILES
system("rm tRNAstructurenumber.txt"); 
open(out1, ">tRNAstructurenumber.txt") or die("Cannot open output file: tRNAref_andvariants.fasta, USAGE: tRNAvariantScript.pl GtRNAdb.txt variants.txt");

open(inp0, "$structure") or die("Cannot open database file");
my @structureline;
my $structurei = 0;
my @reference;
my $genei = 0;
my @sequencelines;

while(<inp0>){
	my @inpfield = split " ";
	my $line = $_; #collects everything in one line and stores it in the variable $line
	my $FirstCharTest = substr($inpfield[0],0,1); #Tests if the line starts with a >
	
	if ($FirstCharTest eq '>' or $FirstCharTest eq '.'){
	my @structureline;
	  $structureline[$structurei] = $line;
	  $structurei++;
	}
	
	else{
	  #Counter keeps track of how many genes
	  my @seq = split("", $inpfield[0]);
	  push @inpfield, \@seq;
	  $reference[$genei] = \@inpfield;
	  #Stores the sequence in a variable for that gene #
	  push @sequencelines, $reference[$genei];
	  $genei++;
	}
}
close(inp0);
my $total = ($genei-1);

open (inp1, "$numbering") or die("Cannot open numbering file");
	my %aminoacidstructure;
	my %NumberofPositions;
	my @aminoacidstructure;
	my $aacount = 0;
	my @aminoacid;
	
while(<inp1>){
	my $line = $_;
	my @inpfield = split(' ', $line);
	$aminoacid[$aacount] = shift @inpfield;
	$aminoacidstructure[$aacount] = \@inpfield;
	$NumberofPositions{$aminoacid[$aacount]} = @inpfield;
	$aacount++;
}

my @numbering;
my $i;
my $counter;


for ($genei = 0; $genei <= $total; $genei++){
	
	for (my $j = 0; $j < $aacount; $j++){	
		
		my @aminoacidcheck = split("-", $sequencelines[$genei][1]);
		my $aminoacidcheck1 = $aminoacidcheck[1];
		
		if ($aminoacid[$j] eq $aminoacidcheck1){
			$counter = 0;
			@numbering = ();
			push @numbering, $sequencelines[$genei][1];
	
			for ($i = 0; $i < $NumberofPositions{"$aminoacidcheck1"}; $i++){
				if ($sequencelines[$genei][4][$i] ne "." && $sequencelines[$genei][4][$i] ne "X"){
					if ($aminoacidstructure[$j][$i] eq "17A" && $numbering[$counter] ne "17"){
					push @numbering, "17";
					$counter++;
					}
					if ($aminoacidstructure[$j][$i] eq "20A" && $numbering[$counter] ne "20"){
					push @numbering, "20";
					$counter++;
					}
					if ($aminoacidstructure[$j][$i] eq "20B" && $numbering[$counter] ne "20A"){
					push @numbering, "20A";
					$counter++;
					}
					else{		
					push @numbering, $aminoacidstructure[$j][$i];
					$counter++;
					}
				}
			}
			
			
			$, = " ";
			print out1 @numbering;
			print out1 "\n";

		}
		else{
		}
	}
}
