#!/usr/bin/perl
use strict;
use warnings;

#By: Matthew Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

## This is a script to extract the tRNAscan Infernal and eufinder tRNA scores and compare reference score to variant scores
## Input 1: GtRNAdb_eufind.out
## Input 2: GtRNAdb_infernal.out
## Input 3: Variants_eufind.out
## Input 4: Variants_infernal.out

# Copyright (C) 2019 Matthew Berg

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# Reads input files
my $RefEufind = $ARGV[0];
my $RefInfernal = $ARGV[1];
my $VarEufind = $ARGV[2];
my $VarInfernal = $ARGV[3];

# Output files
system("rm eufind_scores.txt");
open(out1, ">eufind_scores.txt") or die("Cannot open output file 1");
system("rm infernal_scores.txt");
open(out2, ">infernal_scores.txt") or die("Cannot open output file 2");

# Saves scores for each tRNA in a hash
open(inp0, "$RefEufind") or die("Cannot open reference eufinder tRNA scores");
	my %EufindRef;

while(<inp0>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq 'H'){
		my @splitLine = split("\t", $_);
		my $tRNA = substr($splitLine[0], 13);
		my $score = $splitLine[8];
		
		$EufindRef{$tRNA} = $score;
	}	

}

close(inp0);

open(inp1, "$RefInfernal") or die("Cannot open reference Infernal tRNA scores");
	my %InfernalRef;
	my %pseudoRef;

while(<inp1>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq 'H'){
		my @splitLine = split("\t", $_);
		my $tRNA = substr($splitLine[0], 13);
		$tRNA =~ s/\s+$//;
		my $score = $splitLine[8];
		my $pseudocheck = 0;
		
		if(scalar @splitLine == 9){
			$pseudocheck = 1;
		}
		
		$InfernalRef{$tRNA} = $score;
		$pseudoRef{$tRNA} = $pseudocheck;
		
	}
		
}

close(inp1);

# Saves the scores for each variant tRNA in a hash

open(inp2, "$VarEufind") or die("Cannot open variant eufinder tRNA scores");
	my %EufindVarScore;
	my %EufindVartRNA;

while(<inp2>){
	chomp;
	
	
	my $check = substr($_, 0, 1);
	
	if($check eq 't'){
	
	my @splitLine = split("\t", $_);
	my @temptRNA = split(/\./, $splitLine[0]);
	my $tRNA = shift @temptRNA;
	my $mutation = "@temptRNA";
	my $score = $splitLine[8];
	
	$EufindVarScore{$mutation} = $score;
	$EufindVartRNA{$mutation} = $tRNA;

	}
}

close(inp2);

open(inp3, "$VarInfernal") or die("Cannot open variant infernal tRNA scores");
	my %InfernalVarScore;
	my %InfernalVartRNA;
	my %pseudoVar;

while(<inp3>){
	chomp;
	
	my $check = substr($_, 0, 1);
	
	if($check eq 't'){
	
	
		my @splitLine = split("\t", $_);
		my @temptRNA = split(/\./, $splitLine[0]);
		
		my $tRNA = shift @temptRNA;
		my $mutation = "@temptRNA";
		my $score = $splitLine[8];
		
		my $pseudocheck = 0;
		
		if(scalar @splitLine == 9){
			$pseudocheck = 1;
		}
		
		$InfernalVarScore{$mutation} = $score;
		$InfernalVartRNA{$mutation} = $tRNA;
		$pseudoVar{$mutation} = $pseudocheck;
		
	}
}

close(inp3);

# Organizes reference scores and variant scores into one output file

print out1 "tRNA\tMutation\tRef_Eufinder\tVar_Eufinder\tDelta_Eufinder\n";

foreach my $keys (keys %EufindVartRNA){
	my $tRNA = $EufindVartRNA{$keys};

	if(defined $EufindRef{$tRNA} and defined $EufindVarScore{$keys}){
		print out1 "$tRNA\t$keys\t$EufindRef{$tRNA}\t$EufindVarScore{$keys}\t";
		my $eufindDelta = $EufindVarScore{$keys} - $EufindRef{$tRNA};
		print out1 "$eufindDelta\n";
	}
}

print out2 "tRNA\tMutation\tRef_Infernal\tVar_Infernal\tDelta_Infernal\tPseudo_Change\n";
	
foreach my $keys (keys %InfernalVartRNA){	
	my $tRNA = "$InfernalVartRNA{$keys}";
	
	print out2 "$tRNA\t$keys\t$InfernalRef{$tRNA}\t$InfernalVarScore{$keys}\t";
	my $infernalDelta =  $InfernalVarScore{$keys} - $InfernalRef{$tRNA};
	print out2 "$infernalDelta\t";
	
	my $pseudo = $pseudoVar{$keys} - $pseudoRef{$tRNA};
	my $change = "None";
	
	if($pseudo == 1){
		$change = "Change to NOT PSEUDO";
	}
	
	elsif($pseudo == -1){
		$change = "Change to PSEUDO";
	}
	
	print out2 "$change\n";

}
