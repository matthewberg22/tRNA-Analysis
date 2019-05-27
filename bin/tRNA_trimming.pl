#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

#Script to trim off flanking tRNA sequence, count allele frequency, and output mutant sequence
#To run enter tRNA_trimming.pl tRNA_variant_sequence.fasta
#tRNA_variant_sequence.fasta comes from BLAST_analysis.pl and there is one file for every sample (run with tRNA_trimming_run.sh to analyze and merge all files)

#####Reads in input files
##############################################################################

my $tRNAsequence = $ARGV[0];
my $identifier = $ARGV[1];

#####OUTPUT FILES
##############################################################################
open(out1, ">>total_alleles.txt") or die("Cannot open output1 file");
open(out2, ">>tRNAmutants.txt") or die("Cannot open output2 file");
open(out3, ">>problem_loci.txt") or die("Cannot open output3 file");  

#####Reads in tRNA variant sequences
###############################################################################

print out1 "# $identifier \n";

open(inp0, "$tRNAsequence") or die("Cannot open tRNA_sequence file");
	my $tRNAname;
	my @tRNA;
	my @status;
	my @mutation;
	my @coverage;
	my @sequence;
	my @genomeposition;
	my $sequenceCounter = 0;
	my %Ungapped;
	my %Location;

while(<inp0>){
	chomp;
	my $check = substr($_, 0, 1);
	
	if($check eq ">"){
		my @splitLine = split("\t", $_);
		$tRNAname = substr($splitLine[0], 1);
		$tRNA[$sequenceCounter] = $tRNAname;
		$status[$sequenceCounter] = $splitLine[1];
		my @Tempcoverage = split(":", $splitLine[2]);
		$coverage[$sequenceCounter] = $Tempcoverage[1];		

		
		if($status[$sequenceCounter] eq 'MUT'){
			my @mutationlist;
			$genomeposition[$sequenceCounter] = $splitLine[3];
			
			my @Tempstrand = split('\(', $genomeposition[$sequenceCounter]);
			my $strand = $Tempstrand[1];
			my @Tempstart = split(":", $Tempstrand[0]);
			my @Tempstart2 = split("-", $Tempstart[1]);
			my $start = $Tempstart2[0];
			my $end = $Tempstart2[1];
			
			my @mutationholder = split(" ", $splitLine[4]);
			foreach(@mutationholder){
				my @splitMut = split(":", $_);
				
				if($tRNAname eq 'tRNA-Glu-TTC-3-1'){
				
					if($splitMut[1] < $start + 25 or $splitMut[1] > $end - 5){
					}
					
					else{
					push @mutationlist, $_;
					}
				}
				
				elsif($tRNAname eq 'tRNA-Glu-TTC-4-1'){
				
					if($splitMut[1] > $end - 25 or $splitMut[1] < $start + 5){
					}
					
					else{
					push @mutationlist, $_;
					}
				}
				else{
					if(($splitMut[1] < $start + 20 and $strand eq '+)') or ($splitMut[1] > $end - 5 and $strand eq '+)')){
					}
					elsif(($strand eq '-)' and $splitMut[1] > $end - 20) or ($splitMut[1] < $start + 5 and $strand eq '-)')){
					}
					else{
						push @mutationlist, $_;
					}
				}
			}
			
			my $changeCheck = scalar @mutationlist;
			
			if($changeCheck == 0){
				$status[$sequenceCounter] = 'WT';
			}
			
			elsif($changeCheck >= 1){
				$mutation[$sequenceCounter] = join(" ", @mutationlist);
			}
		}
	}	
	
	else{
	
		if($tRNAname eq 'tRNA-Glu-TTC-3-1' or $tRNAname eq 'tRNA-Glu-TTC-4-1'){
			$sequence[$sequenceCounter] = substr($_, 25);
			$sequence[$sequenceCounter] = substr($sequence[$sequenceCounter], 0, -5);
			$sequenceCounter++;
		
			my $GapTest = substr($_, 25);
			$GapTest = substr($GapTest, 0, -5);
		
			if(index($GapTest, "-") != -1){
				$GapTest =~ tr/-//d;
				$Ungapped{$GapTest}++;
				$Location{$sequenceCounter} = $GapTest;
			}
		}
		
		else{
			$sequence[$sequenceCounter] = substr($_, 20);
			$sequence[$sequenceCounter] = substr($sequence[$sequenceCounter], 0, -5);
			$sequenceCounter++;
		
			my $GapTest = substr($_, 20);
			$GapTest = substr($GapTest, 0, -5);
		
			if(index($GapTest, "-") != -1){
			$GapTest =~ tr/-//d;
			$Ungapped{$GapTest}++;
			$Location{$sequenceCounter} = $GapTest;
			}
		}
	}
				
}

foreach my $keys (keys %Ungapped){
	if ($Ungapped{$keys} > 1){
	my @matching_keys = grep { $Location{$_} eq $keys} keys %Location;
	my $counter = 1;
	my $keepPosition;
	foreach(@matching_keys){
		if($counter == 1){
		$keepPosition = $_;
		$counter++;
		}
		elsif($counter > 1){
		splice @tRNA, $_, 1;
		$coverage[$_] = $coverage[$_] + $coverage[$keepPosition];
		splice @sequence, $_, 1;
		splice @status, $_, 1;
		splice @mutation, $_, 1; 
		}
	}
	}
}

my $counter = 0;
my $mutcounter = 0;
my %tRNAs;
my %seq;
my %numberAlleles;

#Hashes them  down after they are trimmed (ie if there is a mutation in the 20 bp flanking, but the rest is WT it merges that with the WT sequence)
foreach(@tRNA){
	my $name;
	if($status[$counter] eq 'WT'){
	$name = 'WT';	
	}
	elsif($status[$counter] eq 'MUT'){
	$name = $mutation[$counter];
	}
	
	$tRNAs{$_}{$name} = 0;
	$seq{$_}{$name} = $sequence[$counter];
	$counter++;
}

$counter = 0;

#Re-counts the coverage after they have been hashed down
foreach(@tRNA){
	my $name;
	if($status[$counter] eq 'WT'){
	$name = 'WT';	
	}
	elsif($status[$counter] eq 'MUT'){
	$name = $mutation[$counter];
	
	}

	$tRNAs{$_}{$name} = $tRNAs{$_}{$name} + $coverage[$counter];
	$counter++;
}

#Determines if there is one copy of the allele or two (assumes if we only see one allele at a loci it is present in two copies)
my %copynumber;

foreach my $tRNA (sort keys %tRNAs){
	if((scalar %{$tRNAs{$tRNA}}) == 1){
		foreach my $names (sort keys %{$tRNAs{$tRNA}}){
		$copynumber{$tRNA}{$names} = 2;
		}
	}
	
	else{
		foreach my $names (sort keys %{$tRNAs{$tRNA}}){
		$copynumber{$tRNA}{$names} = 1;
		}
	}
}


#Counts the number of alleles
foreach my $tRNA (sort keys %tRNAs){
	foreach my $names (sort keys %{$tRNAs{$tRNA}}){
	
	$numberAlleles{$tRNA}++;
	
	#Counts the number of tRNAs that deviate from the WT
	if($names ne 'WT'){
	$mutcounter++;
	
	print out2 "# $identifier\n";
	print out2 ">$tRNA\t";
	print out2 "$names\tCOV:$tRNAs{$tRNA}{$names}\t$copynumber{$tRNA}{$names}\n$seq{$tRNA}{$names}\n\n";
	}
	
	}
}

#Loops through all of the allele counts and annotates them assuming a diploid (humans)
my %allele_output;


foreach my $key (sort keys %numberAlleles){
	$allele_output{$key} = 0;
	my @tRNALength = split(" ", $key);
	
	if(scalar(@tRNALength) > 1){
		$allele_output{$key} = 2 * scalar(@tRNALength);
	}
	
	elsif(scalar(@tRNALength) == 1){
		if($numberAlleles{$key} > 2){
			$allele_output{$key} = $numberAlleles{$key};
			
			print out3 "$identifier\n";
			print out3 "$key\t$numberAlleles{$key}\n";
			foreach my $names (sort keys %{$seq{$key}}){
				print out3 "$names\tCOV:$tRNAs{$key}{$names}\n $seq{$key}{$names}\n";
			}
			print out3 "\n";
		}
		
		else{
			$allele_output{$key} = 2;
		}
	}
}

foreach my $keys (sort keys %allele_output){
	print out1 "$keys\t$allele_output{$keys}\n";
}

