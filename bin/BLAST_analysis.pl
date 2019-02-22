#!/usr/bin/perl 
use strict;
use warnings;

use Scalar::Util qw(looks_like_number);

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#January 22, 2019

#Script to parse through BLAST output, select hits that contain full length tRNAs plus flanking sequence and sort out possibly ambiguous hits
#To run enter script.pl BLAST.blast IDENTIFIER_Number Coverage_Cutoff_Number
#BLAST.blast contains the output of a BLAST with reference tRNA as query and all the paired end reads as the database
#IDENTIFIER_Number is the unique identifier for the sample you are analyzing
#Coverage_Cutoff_Number is the cutoff you want to use (if you don't know the cutoff and want to get everything set this as 0)

#####Reads in input files
##############################################################################

my $BLAST = $ARGV[0]; 
my $PatientNo = $ARGV[1];
my $coveragecutoff = $ARGV[2];

######OUTPUT FILES
##############################################################################

system("rm BLAST_analysis_$PatientNo.txt"); 
open(out1, ">BLAST_analysis_$PatientNo.txt") or die("Cannot open output1 file");
system("rm tRNA_sequence-$PatientNo-Cov$coveragecutoff.fasta");
open(out2, ">tRNA_sequence-$PatientNo-Cov$coveragecutoff.fasta") or die("Cannot open output2 file");
system("rm Allcoverage_$PatientNo.txt");
open(out3, ">Allcoverage_$PatientNo.txt") or die("Cannot open output2 file");

#####Reads in BLAST file and saves the information
##############################################################################

open(inp0, "$BLAST") or die("Cannot open BLAST file");
              my @tRNAUpCoord; #The genomic 5' co-ordinate of the tRNA with flanking sequence
              my @tRNADnCoord; #The genomic 3' co-ordinate of the tRNA with flanking sequence
              my @tRNARefname; #The name of the reference tRNA from the genome
              my @tRNAchr; #The chromosome the reference tRNA is on
              my @strand;
	      my %StartStop; #A hash to store the start and stop coordinates of the hit within the read
	      my %allele; #A hash to store the btop output
	      my %seq; #A hash to store the sequence of the BLAST hit (read)
	      my %tRNAqueryLength; #A hash to store the total length of each reference tRNA (including 20 bp flanking 5')

#Reads line by line through the BLAST output file              
while(<inp0>){
    chomp;

    my $check = substr($_,0,1); #Tests if the line starts with # meaning it is not a BLAST hit
              
    #If the line is an information line, check if the line has information about the Query (starts with # Query) and store tRNA name and co-ordinates
    if($check eq "#"){
        my @splitLine = split(" ", $_); #Splits the line on spaces
	my $check2 = substr($splitLine[1],0,1);
        
	#If the first word after the # starts with Q (ie its Query) stores the tRNA name and co-ordinates
	if($check2 eq "Q"){
            push @tRNARefname, substr($splitLine[2],13);
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
		my $tRNAname = substr($splitLine[0],13);
		$tRNAqueryLength{$tRNAname} = $splitLine[14];
		push @{ $StartStop{$splitLine[1]}{$tRNAname} }, $splitLine[8].",".$splitLine[9];
		push @{ $allele{$splitLine[1]}{$tRNAname} }, $splitLine[12];
		push @{ $seq{$splitLine[1]}{$tRNAname} }, $splitLine[13];
        }
    }
}

#Creates four hashs, storing information about the reference tRNAs and where they are in the genome
my %tRNAstrand;
@tRNAstrand{@tRNARefname} = @strand;
my %UpCoord;
@UpCoord{@tRNARefname} = @tRNAUpCoord;
my %DnCoord;
@DnCoord{@tRNARefname} = @tRNADnCoord;
my %Chr;
@Chr{@tRNARefname} = @tRNAchr;

close(inp0); 

#####Counts how many reads contain a full length tRNA and how many reads contain two 'duplicated'tRNAs
#####################################################################################################
my $totalfulllength = 0;
my $twotRNAreads = 0;

foreach my $readID (keys %StartStop){
	$totalfulllength++;
	foreach my $tRNAz (keys %{$StartStop{$readID}}){
		if( $#{$StartStop{$readID}{$tRNAz}} >=1 ){
			$twotRNAreads++;
		}
	}
}

#####Checks if readID are unique or are hits for multiple tRNAs
##############################################################################

my $MultipleMappingReadsCounter = 0;
my %confidentalleles;
my $startstop;
my $allele;
my $seq;
my $unresolvable = 0;
my $resolvable = 0;
my %resolvabletRNAs;

####Checks if any of the reads have two tRNAs in them and if true removes the readID from the list
############################################################################################################

#Loops through all the keys of alleles (ie all the readIDs that contain a full length tRNA)
foreach my $readID (sort keys %allele){
	
	my $trip = 0; #Counter to keep track of which readIDs contain more than one full length tRNA hit
	
	#For each tRNA that is a hit for the readID, are there >1 start and stop locations for the same tRNA? If yes, trip counter is set to 1
	foreach my $tRNAz (keys %{$StartStop{$readID}}){
		if( $#{$StartStop{$readID}{$tRNAz}} >=1 ){
			$trip = 1;
		}
	}
	
	#If trip counter is 1, that means a read has a 'duplicated' tRNA, which is an artifact of merging the paired end reads. Therefore the read is deleted
	if($trip == 1){

	delete $StartStop{$readID};
	delete $seq{$readID};
	delete $allele{$readID};
	
	}
}

####Determines which reads are hits for only one reference tRNA and which are hits for multiple tRNAs and need to be looked at closer to be resolved
###############################################################################################################

my $UniqueReadsCounter = 0;
my $check_tRNA;


foreach my $readID (sort keys %allele){
	my %allelecheck;
	my @tRNAfamily;

	my $mappings = scalar keys %{ $allele{$readID} };
	
	#If there is only one tRNA associated with each readID, it is a unique mapping and doesn't need to be resolved further
	if($mappings == 1){
		$UniqueReadsCounter++;
		for my $key (keys %{$allele{$readID}}){
		$check_tRNA = $key;
		}
		
		if ($check_tRNA eq 'tRNA-Phe-GAA-10-1'){

			my $check_seq = $seq{$readID}{$check_tRNA}[0];
			my $check_nuc = substr($check_seq, 53, 1);
			
			
			if($check_nuc eq 'A'){
				delete $StartStop{$readID};
				delete $seq{$readID};
				delete $allele{$readID};
			}
			
		}
				
	}
	
	#If the read is associated with mutltiple tRNAs, it needs to be sorted out
	elsif($mappings > 1){
		$MultipleMappingReadsCounter++;
		
		#Make a list of all the tRNAs that could be the read
		foreach my $tRNAz (keys %{$allele{$readID}}){
			push @tRNAfamily, $tRNAz;
					
			$allelecheck{$allele{$readID}{$tRNAz}[0]} = 0; #A hash where the keys are all the btop outputs for each tRNA corresponding to the read

			$startstop = $StartStop{$readID}{$tRNAz}[0];
			$seq = $seq{$readID}{$tRNAz}[0];
			$allele = $allele{$readID}{$tRNAz}[0];
		
		}
		
		#If the alleles associated with each tRNA for a single read are all the same, we cannot resolve which reference tRNA the read corresponds to
		if(scalar keys %allelecheck == 1){
			$unresolvable++; #Keeps track of how many tRNA groups cannot be resolved
			my $mergedtRNA = join(" ", sort(@tRNAfamily)); #Merge all the tRNAs that it could be into one string
			delete $StartStop{$readID}; #Delete all the information for this readID
			delete $seq{$readID};
			delete $allele{$readID};
			push @{$StartStop{$readID}{$mergedtRNA}}, $startstop; #Rewrite the information with this readID with the merged tRNA as the tRNA name
			push @{$seq{$readID}{$mergedtRNA}}, $seq;
			push @{$allele{$readID}{$mergedtRNA}}, $allele;
		}
		
		#If the alleles associated with each tRNA for a single read are different, we might be able to tell which reference tRNA the read truly corresponds to
		elsif(scalar keys %allelecheck > 1){
			my $mergedtRNA = '';
			my $tripcounter = 0;
			my @WTtRNA;
			my $keytRNA;
			
			#For each tRNA in the family, if one is WT it is probably the tRNA that truly maps to the read
			foreach(sort @tRNAfamily){
				if($tRNAqueryLength{$_} eq $allele{$readID}{$_}[0]){
					push @WTtRNA, $_;
					$keytRNA = $_;
					$tripcounter = 1;
				}
			}
			
			if($tripcounter == 1){
				$mergedtRNA = join(" ", sort(@WTtRNA));
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
			
			#If no tRNA corresponds to a WT match, then using differences that we have determined by hand between similar tRNAs, we might be able to tell which tRNA corresponds to the read
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
					my $check_nuc = substr($check_seq, 0, 10);	
					
					if($check_nuc eq 'AATAAAGCTC'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-1";
					}
					elsif($check_nuc eq 'AATAATGCTG'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-10-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-10-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-10-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-10-1";
						}
					elsif($check_nuc eq 'AATAATGCTC'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-3";
					}
										
				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-13-2/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-12-1/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-12-3/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
					my $check_nuc = substr($check_seq, 0, 19);
					
					if($check_nuc eq 'AATAAAGCTCCCCTGATGT'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-1"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-1";
					}
					elsif($check_nuc eq 'AATAATGCTCCCCTGATGT'){
						$startstop = $StartStop{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$seq = $seq{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-12-3"}[0];
						$mergedtRNA = "tRNA-Ala-AGC-12-3";
					}
					elsif($check_nuc eq 'AATAAAGCTCCCCTGACGT'){
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
						$seq = $seq{$readID}{"tRNA-Ala-AGC-16-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ala-AGC-16-1"}[0];
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
				
				
				
				elsif($tRNAfamilystring =~ m/tRNA-Ala-AGC-12-2/ && $tRNAfamilystring =~ m/tRNA-Ala-AGC-14-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ala-AGC-12-2"}[0];
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
						$allele = $allele{$readID}{"tRNA-Gly-CCC-6-1"}[0];
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
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-10-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-12-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-10-1"}[0];
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

				
				
				elsif($tRNAfamilystring =~ m/tRNA-Asn-GTT-10-1/ && $tRNAfamilystring =~ m/tRNA-Asn-GTT-8-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Asn-GTT-10-1"}[0];
					my $check_nuc = substr($check_seq, 72, 1);
					
					
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
					my $check_nuc = substr($check_seq, 41, 1);					
					
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
				
				elsif($tRNAfamilystring =~ m/nmt-tRNA-Gln-TTG-7-1/ && $tRNAfamilystring =~ m/nmt-tRNA-Gln-TTG-13-1/){
					
					my $check_seq = $seq{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
					my $check_nuc = substr($check_seq, 18, 1);
					
					
					if($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"nmt-tRNA-Gln-TTG-13-1"}[0];
						$seq = $seq{$readID}{"nmt-tRNA-Gln-TTG-13-1"}[0];
						$allele = $allele{$readID}{"nmt-tRNA-Gln-TTG-13-1"}[0];
						$mergedtRNA = "nmt-tRNA-Gln-TTG-13-1";
					}
					elsif($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
						$seq = $seq{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
						$allele = $allele{$readID}{"nmt-tRNA-Gln-TTG-7-1"}[0];
						$mergedtRNA = "nmt-tRNA-Gln-TTG-7-1";
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
				
				elsif($tRNAfamilystring =~ m/tRNA-Met-CAT-5-1/ && $tRNAfamilystring =~ m/tRNA-Met-CAT-7-1/){
					
					my $check_seq = $seq{$readID}{"tRNA-Met-CAT-5-1"}[0];
					my $check_nuc = substr($check_seq, 6, 1);
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Met-CAT-5-1"}[0];
						$seq = $seq{$readID}{"tRNA-Met-CAT-5-1"}[0];
						$allele = $allele{$readID}{"tRNA-Met-CAT-5-1"}[0];
						$mergedtRNA = "tRNA-Met-CAT-5-1";
					}
					elsif($check_nuc eq 'A'){
						$startstop = $StartStop{$readID}{"tRNA-Met-CAT-7-1"}[0];
						$seq = $seq{$readID}{"tRNA-Met-CAT-7-1"}[0];
						$allele = $allele{$readID}{"tRNA-Met-CAT-7-1"}[0];
						$mergedtRNA = "tRNA-Met-CAT-7-1";
					}

				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Ile-AAT-1-1/ && $tRNAfamilystring =~ m/tRNA-Ile-AAT-7-1/ && $tRNAfamilystring =~ m/tRNA-Ile-AAT-7-2/){
					
					my $check_seq = $seq{$readID}{"tRNA-Ile-AAT-1-1"}[0];
					my $check_nuc = substr($check_seq, 3, 1);
					
					if($check_nuc eq 'T'){
						$startstop = $StartStop{$readID}{"tRNA-Ile-AAT-1-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ile-AAT-1-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ile-AAT-1-1"}[0];
						$mergedtRNA = "tRNA-Ile-AAT-1-1";
					}
					elsif($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Ile-AAT-7-1"}[0];
						$seq = $seq{$readID}{"tRNA-Ile-AAT-7-1"}[0];
						$allele = $allele{$readID}{"tRNA-Ile-AAT-7-1"}[0];
						$mergedtRNA = "tRNA-Ile-AAT-7-1 tRNA-Ile-AAT-7-2";
					}

				}
				
				elsif($tRNAfamilystring =~ m/tRNA-Leu-AAG-1-1/ && $tRNAfamilystring =~ m/tRNA-Leu-AAG-1-2/ && $tRNAfamilystring =~ m/tRNA-Leu-AAG-1-3/){
					
					my $check_seq = $seq{$readID}{"tRNA-Leu-AAG-1-1"}[0];
					my $check_nuc = substr($check_seq, 8, 1);
					
					if($check_nuc eq 'G'){
						$startstop = $StartStop{$readID}{"tRNA-Leu-AAG-1-1"}[0];
						$seq = $seq{$readID}{"tRNA-Leu-AAG-1-1"}[0];
						$allele = $allele{$readID}{"tRNA-Leu-AAG-1-1"}[0];
						$mergedtRNA = "tRNA-Leu-AAG-1-1 tRNA-Leu-AAG-1-2";
					}
					elsif($check_nuc eq 'C'){
						$startstop = $StartStop{$readID}{"tRNA-Leu-AAG-1-3"}[0];
						$seq = $seq{$readID}{"tRNA-Leu-AAG-1-3"}[0];
						$allele = $allele{$readID}{"tRNA-Leu-AAG-1-3"}[0];
						$mergedtRNA = "tRNA-Leu-AAG-1-3";
					}

				}
			
				else{;
				}
				
				if($mergedtRNA eq ''){
					delete $StartStop{$readID};
					delete $seq{$readID};
					delete $allele{$readID};
				}
				
				elsif($mergedtRNA ne ''){
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
}


#####Collects stats on how many tRNAs were identified, coverage for each allele, etc.
##############################################################################

my %counts; #A hash to correlate tRNA-allele (key) and number of times its seen/coverage (value)
my %totaltRNAcoverage; #A hash to store each tRNA (key) and how many times reads were seen that map to that tRNA (value)
my %totaltRNAalleleCounter;
my %StoreSeq; #A hash to store each tRNA-allele (key) and a sequence (value)

foreach(@tRNARefname){
	$totaltRNAalleleCounter{$_} = 0;
}
			
foreach my $readID (keys %allele){
	foreach my $tRNAz (keys %{$allele{$readID}}){
	
	my $uniqueallele = $allele{$readID}{$tRNAz}[0];	
	$counts{$tRNAz}{$uniqueallele}++;
	$totaltRNAcoverage{$tRNAz}++;
	$StoreSeq{$tRNAz}{$uniqueallele} = $seq{$readID}{$tRNAz}[0];
	
	
	}
}


my $uniquesequences = 0; #A counter for the total number of unique sequences with >10x coverage seen

#Filters out anything with coverage less than 10x
foreach my $tRNAz (keys %counts){
	for my $uniqueallele (keys %{$counts{$tRNAz}}){
	
	print out3 "$PatientNo\t$tRNAz\t$uniqueallele\t$counts{$tRNAz}{$uniqueallele}\n";
	
		if($counts{$tRNAz}{$uniqueallele} > $coveragecutoff){
	
			$uniquesequences++;
			my @brokendown = split(" ", $tRNAz);
			foreach(@brokendown){
				$totaltRNAalleleCounter{$_}++;
			}
	
		}	
	
	}
}

#Determines which tRNAs have been identified and which are missing (with 10x coverage or greater)
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

#Determines how many tRNAs in total have been identified, even those with low coverage (ie less than 10x)
my $totallowcovtRNAs = 0;

foreach my $keys (keys %totaltRNAcoverage){
	my @splitgrouptRNA = split(" ", $keys);
	foreach(@splitgrouptRNA){
	$totallowcovtRNAs++;
	}
}

print out1 "***** IDENTIFER No.: $PatientNo *****\n";
print out1 "Coverage Cutoff: $coveragecutoff\n";
print out1 "Total reads with full length tRNAs: $totalfulllength\n";
print out1 "Total confident reads: $UniqueReadsCounter\n";
print out1 "Total ambiguous reads: $MultipleMappingReadsCounter\n";
print out1 "Reads with two tRNAs: $twotRNAreads\n";
print out1 "Total tRNAs identified: $totaltRNAs ($totallowcovtRNAs)\n";
print out1 "Missing/low coverage tRNAs: $missingtRNAs\n";
print out1 "Total unique sequences: $uniquesequences\n";
print out1 "#############################################################################\n";


####Creates information for second output containing tRNA name, information, WT/MUT and the mutation position, 20bp flanking + tRNA sequence
###########################################################################################################################

my @sepallele;

foreach my $tRNAz (sort keys %counts){
	for my $tops (sort keys %{$counts{$tRNAz}}){
		if($counts{$tRNAz}{$tops} > $coveragecutoff){
			
			my @brokendown = split(" ", $tRNAz);
			
			
			if($tRNAqueryLength{$brokendown[0]} eq $tops){
				if(scalar @brokendown > 1){
					print out2 ">$tRNAz\tWT\tCov:$counts{$tRNAz}{$tops}\n";
					print out2 "$StoreSeq{$tRNAz}{$tops}\n";
				}
				elsif(scalar @brokendown == 1){			
					print out2 ">$tRNAz\tWT\tCov:$counts{$tRNAz}{$tops}\t$Chr{$tRNAz}:$UpCoord{$tRNAz}-$DnCoord{$tRNAz}($tRNAstrand{$tRNAz})\n";
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
							$mutation[$j] = "$mutChr:$startposition:$ref/$alt";
							$startposition++;
							$j++;
							}
						}
						else{
						$ref = substr($_, 0, 1);
						$alt = substr($_, 1, 1);
						$mutation[$j] = "$mutChr:$startposition:$ref/$alt";
						$startposition++;
						$j++;
						}
					}
				
					my $newj2 = 0;
				
					print out2 ">$tRNAz\tMUT\tCov:$counts{$tRNAz}{$BLASTallele}\t$Chr{$usefultRNA}:$UpCoord{$usefultRNA}-$DnCoord{$usefultRNA}($tRNAstrand{$usefultRNA})\t";
					foreach(@mutation){
						print out2 "$mutation[$newj2] ";
						$newj2++;
					}
					print out2 "\n$StoreSeq{$tRNAz}{$BLASTallele}\n";
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
							$mutation[$j] = "$mutChr:$startposition:$ref/$alt";
							$startposition--;
							$j++;
							}
						}
						else{
						$ref = substr($_, 0, 1);
						$ref =~ tr/-ATGC/-TACG/;
						$alt = substr($_, 1, 1);
						$alt =~ tr/-ATGC/-TACG/;
						$mutation[$j] = "$mutChr:$startposition:$ref/$alt";
						$startposition--;
						$j++;
						}
					}
				
					my $newj2 = 0;
				
					print out2 ">$tRNAz\tMUT\tCov:$counts{$tRNAz}{$BLASTallele}\t$Chr{$usefultRNA}:$UpCoord{$usefultRNA}-$DnCoord{$usefultRNA}($tRNAstrand{$usefultRNA})\t";
					foreach(@mutation){
						print out2 "$mutation[$newj2] ";
						$newj2++;
					}
					print out2 "\t";
					print out2 "\n$StoreSeq{$tRNAz}{$BLASTallele}\n";
				}
			}	
		}
	}
}

print out1 "tRNAname\tTotalCoverage\n";

foreach my $keys (sort keys %totaltRNAcoverage){
	print out1 "$keys\t$totaltRNAcoverage{$keys}\n";
}

