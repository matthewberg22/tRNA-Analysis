#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#February 13, 2019

#Script to count how mutations are present at each position of the tRNA
#To run, enter: mutations_tRNAstructure.pl Annotated_nonredundant_mutants.txt
#Annotated_nonredundant_mutants.txt comes from tRNA_variant_annotation.pl script

#####Reads in input files
##############################################################################

my $annotatednonredundant = $ARGV[0];

#####OUTPUT FILES
##############################################################################
system("rm point_one_map.txt"); 
open(out1, ">point_one_map.txt") or die("Cannot open output1 file");
system("rm indel_one_map.txt"); 
open(out2, ">indel_one_map.txt") or die("Cannot open output2 file");
system("rm doublemutation_matrix.txt"); 
open(out3, ">doublemutation_matrix.txt") or die("Cannot open output3 file");

#####Reads in annotated non-redundant tRNA mutations and parses out relavent information
###############################################################################

open(inp0, "$annotatednonredundant") or die("Cannot open annotated nonredundant variant file");
	my %tRNA;
	my %AF;
	my %NoMutations;
	my %Position;
	my %Location;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check ne "#"){
		my @splitLine = split("\t", $_);
		
		my $mutation = $splitLine[1];
		
		$tRNA{$mutation} = $splitLine[0];
		$AF{$mutation} = $splitLine[2];
		$NoMutations{$mutation} = $splitLine[4];
		$Position{$mutation} = $splitLine[5];
		$Location{$mutation} = $splitLine[6];	
		
	}
}

close(inp0);

my %PositionCount1;
my %PositionCount1uncommon;
my %PositionCount1INDEL;
my %PositionCount1INDELuncommon;
my %Mut1;
my %Mut2;

foreach my $keys (keys %NoMutations){

	if($NoMutations{$keys} == 1){
		my @tempMut = split(":", $keys);
		my @tempMut1 = split("/", $tempMut[2]);
		my $ref = $tempMut1[0];
		my $alt = $tempMut1[1];
		
		if($ref ne '-' && $alt ne '-'){
		$PositionCount1{$Position{$keys}}++;
		$PositionCount1uncommon{$Position{$keys}} = 0;
		}
		
		else{
		$PositionCount1INDEL{$Position{$keys}}++;
		$PositionCount1INDELuncommon{$Position{$keys}} = 0;
		}
	}
	
	elsif($NoMutations{$keys} == 2){
		my @tempPosition = split(" ", $Position{$keys});
		$Mut1{$keys} = $tempPosition[0];
		$Mut2{$keys} = $tempPosition[1];

		
	}
}

my $counter =  0;
my $tRNAtotal = 73;

print out3 "N\t";

for($counter; $counter <= $tRNAtotal; $counter++){
	print out3 "$counter\t";
}

print out3 "\n";

for($counter = 0; $counter <= $tRNAtotal; $counter++){
	print out3 "$counter\t";
	for(my $newcounter = 0; $newcounter <= $tRNAtotal; $newcounter++){
		my $checker = 0;
		foreach my $keys (keys %Mut1){
			if($Mut1{$keys} eq $newcounter && $Mut2{$keys} eq $counter){
			$checker = 1;
			}
			
			else{
			}		
		}
		
		if($checker == 1){
			print out3 "1\t";
		}
		else{
			print out3 "0\t";
		}
	
	}
	print out3 "\n";
	

} 

foreach my $keys (keys %NoMutations){
	if($NoMutations{$keys} == 1 && $AF{$keys} < 0.05){
		my @tempMut = split(":", $keys);
		my @tempMut1 = split("/", $tempMut[2]);
		my $ref = $tempMut1[0];
		my $alt = $tempMut1[1];
		
		if($ref ne '-' && $alt ne '-'){
			$PositionCount1uncommon{$Position{$keys}}++;
		}
		
		else{
		$PositionCount1INDEL{$Position{$keys}}++;
		$PositionCount1INDELuncommon{$Position{$keys}} = 0;
		}
	}

}

print out1 "Position\tAllCounts\tUncommonCounts\n";

foreach my $keys (sort keys %PositionCount1){
	if(substr($keys, 0, 1) eq 'e' || substr($keys, 0, 1) eq 'i' || substr($keys, 0, 1) eq ''){
	}
	else{
	print out1 "$keys\t$PositionCount1{$keys}\t";
	print out1 "$PositionCount1uncommon{$keys}\n";
	}
}

print out2 "Position\tAllCounts\tUncommonCounts\n";

foreach my $keys (sort keys %PositionCount1INDEL){
	if(substr($keys, 0, 1) eq 'e' || substr($keys, 0, 1) eq 'i' || substr($keys, 0, 1) eq ''){
	}
	else{
	print out2 "$keys\t$PositionCount1INDEL{$keys}\t";
	print out2 "$PositionCount1INDELuncommon{$keys}\n";
	}
}	
	
