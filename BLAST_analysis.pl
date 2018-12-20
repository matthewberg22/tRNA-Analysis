#!/usr/bin/perl 
use strict;
use warnings;

use Scalar::Util qw(looks_like_number);

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#December 20, 2018

#To run enter script.pl BLAST.blast reads.fasta
#BLAST.blast contains the output of a BLAST with reference tRNA as query and all the reads as the database
#read.fasta contains all the merged paired end reads from the sequencing run

#####Reads in input files
##############################################################################

my $BLAST = $ARGV[0]; 
my $PatientNo = $ARGV[1];

print "reading BLAST ouput file $BLAST \n";

######OUTPUT FILES
##############################################################################

system("rm BLAST_analysis_$PatientNo.txt"); #for Windows computers use del and for UNIX use rm -f
open(out1, ">BLAST_analysis_$PatientNo.txt") or die("Cannot open output1 file");
system("rm tRNA_sequence_$PatientNo.fasta"); #for Windows computers use del and for UNIX use rm -f
open(out2, ">tRNA_sequence_$PatientNo.fasta") or die("Cannot open output2 file");

#####Reads in BLAST file and saves the information
##############################################################################

open(inp0, "$BLAST") or die("Cannot open BLAST file");
              my @tRNAUpCoord; #The genomic 5' co-ordinate of the tRNA with flanking sequence
              my @tRNADnCoord; #The genomic 3' co-ordinate of the tRNA with flanking sequence
              my @tRNARefname; #The name of the reference tRNA from the genome
              my @tRNAchr; #The chromosome the reference tRNA is on
              my @strand;
			  
			  my %StartStop;
			  my %allele;
			  my %seq;
			  my %tRNAqueryLength;
              
while(<inp0>){
    chomp;
	my @inpfield = split ; #Splits on each line

    my $check = substr($inpfield[0],0,1); #Tests if the line starts with # meaning it is not a BLAST hit
              
    #If the line is an information line, check if the line has information about the Query (starts with # Query) and store tRNA name and co-ordinates
    if($check eq "#"){
        my $check2 = substr($inpfield[1],0,1);
        if($check2 eq "Q"){
            my @splitLine = split(" ", $_);
            push @tRNARefname, substr($splitLine[2],11);
			my @tRNACoords = split(":", $splitLine[3]);
            my @tRNACoords3 = split("=", $tRNACoords[0]);
            push @tRNAchr, $tRNACoords3[1];
            my @tRNACoords2 = split("-", $tRNACoords[1]);
            push @tRNAUpCoord, $tRNACoords2[0];
            push @tRNADnCoord, $tRNACoords2[1];
            my @TempStrand = split("=", $splitLine[6]);
            push @strand, $TempStrand[1];                        
        }
    }
	
              
    #If the line is not an information line, check if BLAST hit covers the full tRNA sequence. If its a full length hit, store information about the BLAST hit including readID, alignment length, read start and stop locations
    elsif($check ne "#"){
        my @splitLine = split("\t", $_);
        
		if($splitLine[3] == $splitLine[14]){
			my $tRNAname = substr($splitLine[0],11);
			$tRNAqueryLength{$tRNAname} = $splitLine[14];
			push @{ $StartStop{$splitLine[1]}{$tRNAname} }, $splitLine[8].",".$splitLine[9];
			push @{ $allele{$splitLine[1]}{$tRNAname} }, $splitLine[12];
			push @{ $seq{$splitLine[1]}{$tRNAname} }, $splitLine[13];
        }
    }
}

my %tRNAstrand;
@tRNAstrand{@tRNARefname} = @strand;

my %UpCoord;
@UpCoord{@tRNARefname} = @tRNAUpCoord;

my %DnCoord;
@DnCoord{@tRNARefname} = @tRNADnCoord;

my %Chr;
@Chr{@tRNARefname} = @tRNAchr;


my @duplicatedtRNA;
my %uniqueduptRNA;
my $uniqueCounter = 0;

close(inp0); #Finish reading in BLAST file

my $totalfulllength = 0;
my $twotRNAreads = 0;

foreach my $readID (keys %StartStop){
	$totalfulllength++;
	foreach my $tRNAz (keys %{$StartStop{$readID}}){
		if( $#{$StartStop{$readID}{$tRNAz}} >=1 ){
			push @duplicatedtRNA, $tRNAz;
			$twotRNAreads++;
		}
	}
}



foreach(@duplicatedtRNA){
	$uniqueduptRNA{$_} = 0;
}

my $numberduplicatedtRNAs = scalar keys %uniqueduptRNA;


#####Checks if readID are unique or are hits for multiple tRNAs
##############################################################################

my $funcounter = 0;

my $UniqueReadsCounter = 0;
my $MultipleMappingReadsCounter = 0;
my %confidentalleles;
my $startstop;
my $allele;
my $seq;
my $unresolvable = 0;
my $resolvable = 0;
my %resolvabletRNAs;


foreach my $readID (sort keys %allele){
	
	my $trip = 0;
	
	foreach my $tRNAz (keys %{$StartStop{$readID}}){
		if( $#{$StartStop{$readID}{$tRNAz}} >=1 ){
			$trip = 1;
		}
	}
	
	if($trip == 1){
		foreach my $tRNAz (keys %{$StartStop{$readID}}){
		for(my $i=0; $i <= $#{ $allele{$readID}{$tRNAz} }; $i++){
				my $newreadID = $i.$readID;
				$StartStop{$newreadID}{$tRNAz}[0] = $StartStop{$readID}{$tRNAz}[$i];
				$allele{$newreadID}{$tRNAz}[0] = $allele{$readID}{$tRNAz}[$i];
				$seq{$newreadID}{$tRNAz}[0] = $seq{$readID}{$tRNAz}[$i];
			
			}
		}

	delete $StartStop{$readID};
	delete $seq{$readID};
	delete $allele{$readID};
	
	}
}


foreach my $readID (sort keys %allele){
	my %allelecheck;
	my @tRNAfamily;


	my $mappings = scalar keys %{ $allele{$readID} };
	
	if($mappings == 1){
		$UniqueReadsCounter++;			
	}
	
	elsif($mappings > 1){
		$MultipleMappingReadsCounter++;
		foreach my $tRNAz (keys %{$allele{$readID}}){

			push @tRNAfamily, $tRNAz;
					
			$allelecheck{$allele{$readID}{$tRNAz}[0]} = 0;

			$startstop = $StartStop{$readID}{$tRNAz}[0];
			$seq = $seq{$readID}{$tRNAz}[0];
			$allele = $allele{$readID}{$tRNAz}[0];
		
		}
	
		if(scalar keys %allelecheck == 1){
			$unresolvable++;
			my $mergedtRNA = join(" ", sort(@tRNAfamily));
			delete $StartStop{$readID};
			delete $seq{$readID};
			delete $allele{$readID};
			push @{$StartStop{$readID}{$mergedtRNA}}, $startstop;
			push @{$seq{$readID}{$mergedtRNA}}, $seq;
			push @{$allele{$readID}{$mergedtRNA}}, $allele;
		}

		elsif(scalar keys %allelecheck > 1){
			my $mergedtRNA;
			my $tripcounter = 0;
			my @mertRNA;
			my $keytRNA;
			foreach(sort @tRNAfamily){
				if($tRNAqueryLength{$_} eq $allele{$readID}{$_}[0]){
					push @mertRNA, $_;
					$keytRNA = $_;
					$tripcounter = 1;
				}
			}
			
			if($tripcounter == 1){
				$mergedtRNA = join(" ", @mertRNA);
				$startstop = $StartStop{$readID}{$keytRNA}[0];
				$seq = $seq{$readID}{$keytRNA}[0];
				$allele = $allele{$readID}{$keytRNA}[0];
				delete $StartStop{$readID};
				delete $seq{$readID};
				delete $allele{$readID};
				push @{$StartStop{$readID}{$mergedtRNA}}, $startstop;
				push @{$seq{$readID}{$mergedtRNA}}, $seq;
				push @{$allele{$readID}{$mergedtRNA}}, $allele;
			}
		
			elsif($tripcounter == 0){
				my $tRNAfamilystring = join(" ", sort(@tRNAfamily));
				my $allele = "";
				my $seq = "";
				my $startstop = "";
				
				if($tRNAfamilystring =~ m/tRNA-Ala-AGC-12-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-13-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-13-2/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
					my $check_nuc = substr($check_seq, 3, 1);
								
					if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-1";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-13-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-13-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-13-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-13-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-13-2";
					}
				}
			
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-12-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-12-3/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-10-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
					my $check_nuc = substr($check_seq, 6, 1);
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-1";
					}
					elsif($check_nuc eq 'T'){
						my $check_nuc = substr($check_seq, 10, 1);
						
						if($check_nuc eq 'G'){
							$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-10-1"}[0];
							$seq = $seq{$readID}{"tRNA-Ala-AGC-10-1"}[0];
							$allele = $allele{$readID}{"tRNA-Ala-AGC-10-1"}[0];
							$mergedtRNA = "tRNA-Ala-AGC-10-1";
						}
						elsif($check_nuc eq 'C'){
							$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-3"}[0];
							$seq = $seq{$readID}{"tRNA-Ala-AGC-12-3"}[0];
							$allele = $allele{$readID}{"tRNA-Ala-AGC-12-3"}[0];
							$mergedtRNA = "tRNA-Ala-AGC-12-3";
						}
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-13-2/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-12-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-12-3/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
					my $check_nuc = substr($check_seq, 6, 1);
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-3";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-13-2";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-12-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-12-3/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
					my $check_nuc = substr($check_seq, 6, 1);
					
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-3";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-13-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-13-2/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-13-1"}[0];
					my $check_nuc = substr($check_seq, 3, 1);
					
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-13-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-13-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-13-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-13-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-13-2"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-13-2";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-12-2/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-14-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-16-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-12-2"}[0];
					my $check_nuc = substr($check_seq, 31, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-16-1"}[0];
						my $seq = $seq{$readID}{"tRNA-Ala-AGC-16-1"}[0];
						my $allele = $allele{$readID}{"tRNA-Ala-AGC-16-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-16-1";
					}
					elsif($check_nuc eq 'T'){
						my $check_nuc = substr($check_seq, 76, 1);

						if($check_nuc eq 'G'){
							$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-2"}[0];
							$seq = $seq{$readID}{"tRNA-Ala-AGC-12-2"}[0];
							$allele = $allele{$readID}{"tRNA-Ala-AGC-12-2"}[0];
							$mergedtRNA = "tRNA-Ala-AGC-12-2";
						}
						elsif($check_nuc eq 'A'){
							$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-14-1"}[0];
							$seq = $seq{$readID}{"tRNA-Ala-AGC-14-1"}[0];
							$allele = $allele{$readID}{"tRNA-Ala-AGC-14-1"}[0];
							$mergedtRNA = "tRNA-Ala-AGC-14-1";
						}
						
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-10-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-12-3/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-10-1"}[0];
					my $check_nuc = substr($check_seq, 10, 1);
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-10-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-10-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-10-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-10-1";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-3";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-19-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-21-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-19-1"}[0];
					my $check_nuc = substr($check_seq, 42, 1);
					
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-21-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-21-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-21-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-21-1";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-19-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-19-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-19-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-19-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Leu-CAG-1-1/ && $tRNAfamilystring =~ m/tRNA-Leu-CAG-1-2/ && $tRNAfamilystring =~ m/tRNA-Leu-CAG-1-3/ && $tRNAfamilystring =~ m/tRNA-Leu-CAG-1-4/ && $tRNAfamilystring =~ m/tRNA-Leu-CAG-1-5/){
					
					my $check_seq = $seq{$readID}{"tRNA-Leu-CAG-1-1"}[0];
					my $check_nuc = substr($check_seq, 4, 1);
					
					
					if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Leu-CAG-1-1"}[0];
						$seq = $seq{$readID}{"tRNA-Leu-CAG-1-1"}[0];
						$allele = $allele{$readID}{"tRNA-Leu-CAG-1-1"}[0];
						$mergedtRNA = "tRNA-Leu-CAG-1-1";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Leu-CAG-1-2"}[0];
						$seq = $seq{$readID}{"tRNA-Leu-CAG-1-2"}[0];
						$allele = $allele{$readID}{"tRNA-Leu-CAG-1-2"}[0];
						$mergedtRNA = "tRNA-Leu-CAG-1-2 tRNA-Leu-CAG-1-3 tRNA-Leu-CAG-1-4 tRNA-Leu-CAG-1-5";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Gly-CCC-8-1/ && $tRNAfamilystring =~ m/tRNA-Gly-CCC-chr1-100/){
					
					$startstop = $StartStop{$readID}{"tRNA-Gly-CCC-8-1"}[0];
					$seq = $seq{$readID}{"tRNA-Gly-CCC-8-1"}[0];
					$allele = $allele{$readID}{"tRNA-Gly-CCC-8-1"}[0];
					$mergedtRNA = "tRNA-Gly-CCC-8-1 tRNA-Gly-CCC-chr1-100";
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Gly-CCC-1-1/ && $tRNAfamilystring =~ m/tRNA-Gly-CCC-1-2/){
					
					my $check_seq = $seq{$readID}{"tRNA-Gly-CCC-1-2"}[0];
					my $check_nuc = substr($check_seq, 8, 1);
					
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Gly-CCC-1-1"}[0];
						$seq = $seq{$readID}{"tRNA-Gly-CCC-1-1"}[0];
						$allele = $allele{$readID}{"tRNA-Gly-CCC-1-1"}[0];
						$mergedtRNA = "tRNA-Gly-CCC-1-1";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Gly-CCC-1-2"}[0];
						$seq = $seq{$readID}{"tRNA-Gly-CCC-1-2"}[0];
						$allele = $allele{$readID}{"tRNA-Gly-CCC-1-2"}[0];
						$mergedtRNA = "tRNA-Gly-CCC-1-2";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Gly-TCC-2-2/ && $tRNAfamilystring =~ m/tRNA-Gly-TCC-2-3/ && $tRNAfamilystring =~ m/tRNA-Gly-TCC-2-4/ && $tRNAfamilystring =~ m/tRNA-Gly-TCC-2-5/ && $tRNAfamilystring =~ m/tRNA-Gly-TCC-4-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Gly-TCC-2-2"}[0];
					my $check_nuc = substr($check_seq, 47, 1);
					
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Gly-TCC-4-1"}[0];
						$seq = $seq{$readID}{"tRNA-Gly-TCC-4-1"}[0];
						$allele = $allele{$readID}{"tRNA-Gly-TCC-4-1"}[0];
						$mergedtRNA = "tRNA-Gly-TCC-4-1";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Gly-TCC-2-2"}[0];
						$seq = $seq{$readID}{"tRNA-Gly-TCC-2-2"}[0];
						$allele = $allele{$readID}{"tRNA-Gly-TCC-2-2"}[0];
						$mergedtRNA = "tRNA-Gly-TCC-2-2 tRNA-Gly-TCC-2-3 tRNA-Gly-TCC-2-4 tRNA-Gly-TCC-2-5";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Gly-CCC-6-1/ && $tRNAfamilystring =~ m/tRNA-Val-CAC-8-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Gly-CCC-6-1"}[0];
					my $check_nuc = substr($check_seq, 7, 1);
					
					
					if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Gly-CCC-6-1"}[0];
						$seq = $seq{$readID}{"tRNA-Gly-CCC-6-1"}[0];
						$allele = $allele{$readID}{"tRNA-Gly-CCC-6-11"}[0];
						$mergedtRNA = "tRNA-Gly-CCC-6-1";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Val-CAC-8-1"}[0];
						$seq = $seq{$readID}{"tRNA-Val-CAC-8-1"}[0];
						$allele = $allele{$readID}{"tRNA-Val-CAC-8-1"}[0];
						$mergedtRNA = "tRNA-Val-CAC-8-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-10-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-12-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-8-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-10-1"}[0];
					my $check_nuc = substr($check_seq, 5, 1);
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-8-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-8-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-8-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-8-1";
					}
					elsif($check_nuc eq 'C'){
						my $check_nuc = substr($check_seq, 65, 1);

						if($check_nuc eq 'G'){
							$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-10-1"}[0];
							$seq = $seq{$readID}{"tRNA-Asn-GTT-10-1"}[0];
							$allele = $allele{$readID}{"tRNA-Asn-GTT-10-1"}[0];
							$mergedtRNA = "tRNA-Asn-GTT-10-1";
						}
						elsif($check_nuc eq 'A'){
							$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-12-1"}[0];
							$seq = $seq{$readID}{"tRNA-Asn-GTT-12-1"}[0];
							$allele = $allele{$readID}{"tRNA-Asn-GTT-12-1"}[0];
							$mergedtRNA = "tRNA-Asn-GTT-12-1";
						}
					}

				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-10-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-8-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-10-1"}[0];
					my $check_nuc = substr($check_seq, 5, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-10-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-10-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-10-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-10-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-8-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-8-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-8-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-8-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-15-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-19-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-15-1"}[0];
					my $check_nuc = substr($check_seq, 53, 1);
					
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-15-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-15-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-15-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-15-1";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-19-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-19-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-19-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-19-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-11-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-11-2/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-14-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-11-1"}[0];
					my $check_nuc = substr($check_seq, 76, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-11-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-11-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-11-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-11-1 tRNA-Asn-GTT-11-2";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-14-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-14-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-14-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-14-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-2-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-3-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-3-2/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-2-1"}[0];
					my $check_nuc = substr($check_seq, 10, 1);
					
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-3-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-3-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-3-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-3-1 tRNA-Asn-GTT-3-2";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-2-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-2-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-2-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-2-1";
					}

				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-16-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-16-2/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-16-3/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-16-4/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-16-5/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-16-1"}[0];
					my $check_nuc = substr($check_seq, 17, 1);

					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-16-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-16-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-16-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-16-1 tRNA-Asn-GTT-16-4";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-16-2"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-16-2"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-16-2"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-16-2 tRNA-Asn-GTT-16-3 tRNA-Asn-GTT-16-5";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-13-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-18-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-13-1"}[0];
					my $check_nuc = substr($check_seq, 77, 1);

					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-13-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-13-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-13-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-13-1";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Asn-GTT-18-1"}[0];
						$seq = $seq{$readID}{"tRNA-Asn-GTT-18-1"}[0];
						$allele = $allele{$readID}{"tRNA-Asn-GTT-18-1"}[0];
						$mergedtRNA = "tRNA-Asn-GTT-18-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/nm-tRNA-Tyr-ATA-chr11-9/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr17-33/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr7-6/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr9-10/){
					
					my $check_seq = $seq{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
					my $check_nuc = substr($check_seq, 6, 1);
					
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-ATA-chr11-9";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr17-33"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr17-33"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr17-33"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr17-33";
					}
					elsif($check_nuc eq 'T'){
						my $check_nuc = substr($check_seq, 43, 1);

						if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr7-6";
						}
						elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr9-10"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr9-10"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr9-10"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr9-10";
						}
					}
				}
				
				elsif($tRNAfamilystring =~ m/nm-tRNA-Tyr-ATA-chr11-9/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr17-33/){
					
					my $check_seq = $seq{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
					my $check_nuc = substr($check_seq, 6, 1);
					
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-ATA-chr11-9";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr17-33"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr17-33"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr17-33"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr17-33";
					}
				}
				
				elsif($tRNAfamilystring =~ m/nm-tRNA-Tyr-ATA-chr11-9/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr7-6/){
					
					my $check_seq = $seq{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
					my $check_nuc = substr($check_seq, 0, 1);
					
					
					if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-ATA-chr11-9"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-ATA-chr11-9";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr7-6";
					}
				}
				
				elsif($tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr7-6/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr9-10/){
					
					my $check_seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
					my $check_nuc = substr($check_seq, 19, 1);
					
					
					if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr7-6"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr7-6";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr9-10"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr9-10"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr9-10"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr9-10";
					}
				}
				
				elsif($tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr1-125/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-GTA-chr14-8/){
					
					my $check_seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr1-125"}[0];
					my $check_nuc = substr($check_seq, 2, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr1-125"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr1-125"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr1-125"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr1-125";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-GTA-chr14-8"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-GTA-chr14-8"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-GTA-chr14-8"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-GTA-chr14-8";
					}
				}
				
				elsif($tRNAfamilystring =~ m/nm-tRNA-Tyr-NNN-chr2-20/ && $tRNAfamilystring =~ m/nm-tRNA-Tyr-NNN-chr2-8/){
					
					my $check_seq = $seq{$readID}{"nm-tRNA-Tyr-NNN-chr2-20"}[0];
					my $check_nuc = substr($check_seq, 55, 1);
					
					
					if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-NNN-chr2-20"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-NNN-chr2-20"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-NNN-chr2-20"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-NNN-chr2-20";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"nm-tRNA-Tyr-NNN-chr2-8"}[0];
						$seq = $seq{$readID}{"nm-tRNA-Tyr-NNN-chr2-8"}[0];
						$allele = $allele{$readID}{"nm-tRNA-Tyr-NNN-chr2-8"}[0];
						$mergedtRNA = "nm-tRNA-Tyr-NNN-chr2-8";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Val-CAC-4-1/ && $tRNAfamilystring =~ m/tRNA-Val-CAC-5-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Val-CAC-4-1"}[0];
					my $check_nuc = substr($check_seq, 3, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Val-CAC-4-1"}[0];
						$seq = $seq{$readID}{"tRNA-Val-CAC-4-1"}[0];
						$allele = $allele{$readID}{"tRNA-Val-CAC-4-1"}[0];
						$mergedtRNA = "tRNA-Val-CAC-4-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Val-CAC-5-1"}[0];
						$seq = $seq{$readID}{"tRNA-Val-CAC-5-1"}[0];
						$allele = $allele{$readID}{"tRNA-Val-CAC-5-1"}[0];
						$mergedtRNA = "tRNA-Val-CAC-5-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Val-CAC-11-1/ && $tRNAfamilystring =~ m/tRNA-Val-CAC-11-2/){
					
					my $check_seq = $seq{$readID}{"tRNA-Val-CAC-11-1"}[0];
					my $check_nuc = substr($check_seq, 12, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Val-CAC-11-1"}[0];
						$seq = $seq{$readID}{"tRNA-Val-CAC-11-1"}[0];
						$allele = $allele{$readID}{"tRNA-Val-CAC-11-1"}[0];
						$mergedtRNA = "tRNA-Val-CAC-11-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Val-CAC-11-2"}[0];
						$seq = $seq{$readID}{"tRNA-Val-CAC-11-2"}[0];
						$allele = $allele{$readID}{"tRNA-Val-CAC-11-2"}[0];
						$mergedtRNA = "tRNA-Val-CAC-11-2";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Val-AAC-1-2/ && $tRNAfamilystring =~ m/tRNA-Val-AAC-1-3/){
					
					my $check_seq = $seq{$readID}{"tRNA-Val-AAC-1-2"}[0];
					my $check_nuc = substr($check_seq, 3, 1);
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Val-AAC-1-2"}[0];
						$seq = $seq{$readID}{"tRNA-Val-AAC-1-2"}[0];
						$allele = $allele{$readID}{"tRNA-Val-AAC-1-2"}[0];
						$mergedtRNA = "tRNA-Val-AAC-1-2";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Val-AAC-1-3"}[0];
						$seq = $seq{$readID}{"tRNA-Val-AAC-1-3"}[0];
						$allele = $allele{$readID}{"tRNA-Val-AAC-1-3"}[0];
						$mergedtRNA = "tRNA-Val-AAC-1-3";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Cys-GCA-10-1/ && $tRNAfamilystring =~ m/tRNA-Cys-GCA-18-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Cys-GCA-10-1"}[0];
					my $check_nuc = substr($check_seq, 0, 1);
					
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Cys-GCA-10-1"}[0];
						$seq = $seq{$readID}{"tRNA-Cys-GCA-10-1"}[0];
						$allele = $allele{$readID}{"tRNA-Cys-GCA-10-1"}[0];
						$mergedtRNA = "tRNA-Cys-GCA-10-1";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Cys-GCA-18-1"}[0];
						$seq = $seq{$readID}{"tRNA-Cys-GCA-18-1"}[0];
						$allele = $allele{$readID}{"tRNA-Cys-GCA-18-1"}[0];
						$mergedtRNA = "tRNA-Cys-GCA-18-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ile-GAT-chrX-5/ && $tRNAfamilystring =~ m/tRNA-Ile-GAT-chrX-7/ && $tRNAfamilystring =~ m/tRNA-Ile-GAT-chrX-9/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ile-GAT-chrX-5"}[0];
					my $check_nuc = substr($check_seq, 71, 1);
					
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ile-GAT-chrX-5"}[0];
						$seq = $seq{$readID}{"tRNA-Ile-GAT-chrX-5"}[0];
						$allele = $allele{$readID}{"tRNA-Ile-GAT-chrX-5"}[0];
						$mergedtRNA = "tRNA-Ile-GAT-chrX-5 tRNA-Ile-GAT-chrX-7";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Ile-GAT-chrX-9"}[0];
						$seq = $seq{$readID}{"tRNA-Ile-GAT-chrX-9"}[0];
						$allele = $allele{$readID}{"tRNA-Ile-GAT-chrX-9"}[0];
						$mergedtRNA = "tRNA-Ile-GAT-chrX-9";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ile-GAT-1-1/ && $tRNAfamilystring =~ m/tRNA-Ile-GAT-1-2/ && $tRNAfamilystring =~ m/tRNA-Ile-GAT-1-3/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ile-GAT-1-1"}[0];
					my $check_nuc = substr($check_seq, 16, 1);

					if($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Ile-GAT-1-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ile-GAT-1-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ile-GAT-1-1"}[0];
						$mergedtRNA = "tRNA-Ile-GAT-1-1";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ile-GAT-1-2"}[0];
						$seq = $seq{$readID}{"tRNA-Ile-GAT-1-2"}[0];
						$allele = $allele{$readID}{"tRNA-Ile-GAT-1-2"}[0];
						$mergedtRNA = "tRNA-Ile-GAT-1-2 tRNA-Ile-GAT-1-3";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Gln-CTG-10-1/ && $tRNAfamilystring =~ m/tRNA-Gln-CTG-8-1/ && $tRNAfamilystring =~ m/tRNA-Gln-CTG-8-2/){
					
					my $check_seq = $seq{$readID}{"tRNA-Gln-CTG-10-1"}[0];
					my $check_nuc = substr($check_seq, 42, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Gln-CTG-10-1"}[0];
						$seq = $seq{$readID}{"tRNA-Gln-CTG-10-1"}[0];
						$allele = $allele{$readID}{"tRNA-Gln-CTG-10-1"}[0];
						$mergedtRNA = "tRNA-Gln-CTG-10-1";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Gln-CTG-8-1"}[0];
						$seq = $seq{$readID}{"tRNA-Gln-CTG-8-1"}[0];
						$allele = $allele{$readID}{"tRNA-Gln-CTG-8-1"}[0];
						$mergedtRNA = "tRNA-Gln-CTG-8-1 tRNA-Gln-CTG-8-2";
					}

				}
				
				elsif($tRNAfamilystring =~ m/nmt-tRNA-Gln-TTG-7-1/ && $tRNAfamilystring =~ m/nmt-tRNA-Gln-TTG-9-1/){
					
					my $check_seq = $seq{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
					my $check_nuc = substr($check_seq, 44, 1);
					
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
						$seq = $seq{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
						$allele = $allele{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
						$mergedtRNA = "nmt-tRNA-Gln-TTG-7-1";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"nmt-tRNA-Gln-TTG-9-1"}[0];
						$seq = $seq{$readID}{"nmt-tRNA-Gln-TTG-9-1"}[0];
						$allele = $allele{$readID}{"nmt-tRNA-Gln-TTG-9-1"}[0];
						$mergedtRNA = "nmt-tRNA-Gln-TTG-9-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Glu-TTC-10-1/ && $tRNAfamilystring =~ m/tRNA-Glu-TTC-8-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Glu-TTC-10-1"}[0];
					my $check_nuc = substr($check_seq, 23, 1);
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Glu-TTC-10-1"}[0];
						$seq = $seq{$readID}{"tRNA-Glu-TTC-10-1"}[0];
						$allele = $allele{$readID}{"tRNA-Glu-TTC-10-1"}[0];
						$mergedtRNA = "tRNA-Glu-TTC-10-1";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Glu-TTC-8-1"}[0];
						$seq = $seq{$readID}{"tRNA-Glu-TTC-8-1"}[0];
						$allele = $allele{$readID}{"tRNA-Glu-TTC-8-1"}[0];
						$mergedtRNA = "tRNA-Glu-TTC-8-1";
					}
					else{
						$startstop = $StartStop{$readID}{"tRNA-Glu-TTC-8-1"}[0];
						$seq = $seq{$readID}{"tRNA-Glu-TTC-8-1"}[0];
						$allele = $allele{$readID}{"tRNA-Glu-TTC-8-1"}[0];
						$mergedtRNA = "tRNA-Glu-TTC-8-1 tRNA-Glu-TTC-10-1";
					}
				}
				
				elsif($tRNAfamilystring =~ m/nmt-tRNA-Ser-TGA-1-1/ && $tRNAfamilystring =~ m/nmt-tRNA-Ser-TGA-2-1/){
					
					my $check_seq = $seq{$readID}{"nmt-tRNA-Ser-TGA-1-1"}[0];
					my $check_nuc = substr($check_seq, 30, 1);
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"nmt-tRNA-Ser-TGA-1-1"}[0];
						$seq = $seq{$readID}{"nmt-tRNA-Ser-TGA-1-1"}[0];
						$allele = $allele{$readID}{"nmt-tRNA-Ser-TGA-1-1"}[0];
						$mergedtRNA = "nmt-tRNA-Ser-TGA-1-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"nmt-tRNA-Ser-TGA-2-1"}[0];
						$seq = $seq{$readID}{"nmt-tRNA-Ser-TGA-2-1"}[0];
						$allele = $allele{$readID}{"nmt-tRNA-Ser-TGA-2-1"}[0];
						$mergedtRNA = "nmt-tRNA-Ser-TGA-2-1";
					}

				}
			
				else{
					print $tRNAfamilystring;
					print " - Problem \n";
				}
				
				delete $StartStop{$readID};
				delete $seq{$readID};
				delete $allele{$readID};
				push @{$StartStop{$readID}{$mergedtRNA}}, $startstop;
				push @{$seq{$readID}{$mergedtRNA}}, $seq;
				push @{$allele{$readID}{$mergedtRNA}}, $allele;
		
			}
		
		}
	}
}

my %counts;
my %totaltRNAcoverage;
my %totaltRNAalleleCounter;
my %StoreSeq;

foreach(@tRNARefname){
	$totaltRNAalleleCounter{$_} = 0;
}
			
foreach my $readID (keys %allele){
	foreach my $tRNAz (keys %{$allele{$readID}}){
	$funcounter++;
	
	my $tops = $allele{$readID}{$tRNAz}[0];	
	$counts{$tRNAz}{$tops}++;
	$totaltRNAcoverage{$tRNAz}++;
	$StoreSeq{$tRNAz}{$tops} = $seq{$readID}{$tRNAz}[0];
	
	}
}


my $unsequences = 0;

foreach my $tRNAz (keys %counts){
	for my $tops (keys %{$counts{$tRNAz}}){
		if($counts{$tRNAz}{$tops} > 9){
	
			$unsequences++;
			my @brokendown = split(" ", $tRNAz);
			foreach(@brokendown){
				$totaltRNAalleleCounter{$_}++;
			}
	
		}	
	
	}
}

my $totaltRNAs = 0;
my $missingtRNAs = "";

foreach(@tRNARefname){
	if($totaltRNAalleleCounter{$_} == 0){
	$missingtRNAs = $missingtRNAs." ".$_;
	}
	else{
	$totaltRNAs++;
	}
}

my $totallowcovtRNAs;

foreach my $keys (keys %totaltRNAcoverage){
	my @splitgrouptRNA = split(" ", $keys);
	foreach(@splitgrouptRNA){
	$totallowcovtRNAs++;
	}
}

print out1 "***** IDENTIFER No.: $PatientNo *****\n";
print out1 "Total reads with full length tRNAs: $totalfulllength\n";
print out1 "Total confident reads: $UniqueReadsCounter\n";
print out1 "Total ambiguous reads: $MultipleMappingReadsCounter\n";
print out1 "Reads with two tRNAs: $twotRNAreads\n";
print out1 "No. of duplicated tRNAs: $numberduplicatedtRNAs\n";  
print out1 "List of Duplicated tRNAs: ";
print out1 "\t$_ " for keys %uniqueduptRNA;
print out1 "\nTotal tRNAs identified: $totaltRNAs ($totallowcovtRNAs)\n";
print out1 "Missing/low coverage tRNAs: $missingtRNAs\n";
print out1 "Total unique sequences: $unsequences\n";
print out1 "#############################################################################\n";

my @sepallele;

foreach my $tRNAz (sort keys %counts){
	for my $tops (sort keys %{$counts{$tRNAz}}){
		if($counts{$tRNAz}{$tops} > 9){
			
			my @brokendown = split(" ", $tRNAz);
			
			
			if($tRNAqueryLength{$brokendown[0]} eq $tops){
				if(scalar @brokendown > 1){
					print out2 ">$tRNAz WT Cov:$counts{$tRNAz}{$tops}\n";
					print out2 "$StoreSeq{$tRNAz}{$tops}\n";
				}
				elsif(scalar @brokendown == 1){			
					print out2 ">$tRNAz $Chr{$tRNAz}:$UpCoord{$tRNAz}-$DnCoord{$tRNAz} ($tRNAstrand{$tRNAz}) WT Cov:$counts{$tRNAz}{$tops}\n";
					print out2 "$StoreSeq{$tRNAz}{$tops}\n";
				}
			}
			else{
				my $usefultRNA = $brokendown[0];
				
				my $BLASTallele = $tops;
				$tops =~ s/(\d+)/"$1 "/eg;
				$tops =~ s/(\D+)/"$1 "/eg;
				@sepallele = split(" ", $tops);
				
				
				my $startposition;
				my $mutChr = $Chr{$usefultRNA};
				
				if($tRNAstrand{$usefultRNA} eq '+'){
					$startposition = $UpCoord{$usefultRNA};
					my $j = 0;
					my @mutation;
					
					foreach(@sepallele){
						my $ref;
						my$alt;
						
						if(looks_like_number($_)){
						$startposition = $startposition + $_;
						}
						elsif(length($_) > 2){
							my @splitmultiples = $_ =~ m/../g;
							my $newj = 0;
							foreach(@splitmultiples){
							$ref = substr($_, 0, 1);
							$alt = substr($_, 1, 1);
							$startposition++;
							$mutation[$j] = "$mutChr:$startposition $ref/$alt";
							$j++;
							}
						}
						else{
						$ref = substr($_, 0, 1);
						$alt = substr($_, 1, 1);
						$startposition++;
						$mutation[$j] = "$mutChr:$startposition $ref/$alt";
						$j++;
						}
					}
				
					my $newj2 = 0;
				
					print out2 ">$tRNAz $Chr{$usefultRNA}:$UpCoord{$usefultRNA}-$DnCoord{$usefultRNA} ($tRNAstrand{$usefultRNA}) ";
					foreach(@mutation){
						print out2 "$mutation[$newj2] ";
						$newj2++;
					}
					print out2 "Cov:$counts{$tRNAz}{$BLASTallele}\n$StoreSeq{$tRNAz}{$BLASTallele}\n";
				}
			
			
				elsif($tRNAstrand{$usefultRNA} eq '-'){
					$startposition = $DnCoord{$usefultRNA};
					my $j = 0;
					my @mutation;
					
					foreach(@sepallele){
						my $ref;
						my$alt;
						
						if(looks_like_number($_)){
						$startposition = $startposition - $_;
						}
						elsif(length($_) > 2){
							my @splitmultiples = $_ =~ m/../g;
							my $newj = 0;
							foreach(@splitmultiples){
							$ref = substr($_, 0, 1);
							$ref =~ tr/-ATGC/-TACG/;
							$alt = substr($_, 1, 1);
							$alt =~ tr/-ATGC/-TACG/;
							$startposition--;
							$mutation[$j] = "$mutChr:$startposition $ref/$alt";
							$j++;
							}
						}
						else{
						$ref = substr($_, 0, 1);
						$ref =~ tr/-ATGC/-TACG/;
						$alt = substr($_, 1, 1);
						$alt =~ tr/-ATGC/-TACG/;
						$startposition--;
						$mutation[$j] = "$mutChr:$startposition $ref/$alt";
						$j++;
						}
					}
				
					my $newj2 = 0;
				
					print out2 ">$tRNAz $Chr{$usefultRNA}:$UpCoord{$usefultRNA}-$DnCoord{$usefultRNA} ($tRNAstrand{$usefultRNA}) ";
					foreach(@mutation){
						print out2 "$mutation[$newj2] ";
						$newj2++;
					}
					print out2 "Cov:$counts{$tRNAz}{$BLASTallele}\n$StoreSeq{$tRNAz}{$BLASTallele}\n";
				}
			}	
		}
	}
}

print out1 "#############################################################################\n";
print out1 "tRNAname\tTotalCoverage\n";

foreach my $keys (keys %totaltRNAcoverage){
	print out1 "$keys\t$totaltRNAcoverage{$keys}\n";
}

