#!/usr/bin/perl 
use strict;
use warnings;

#By Matt Berg
#Department of Biochemistry, University of Western Ontario
#December 15, 2018

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

system("del BLAST_analysis_$PatientNo.txt"); #for Windows computers use del and for UNIX use rm -f
open(out1, ">BLAST_analysis_$PatientNo.txt") or die("Cannot open output1 file");
system("del multiply_mapping_tRNAfamiles.txt"); #for Windows computers use del and for UNIX use rm -f
open(out2, ">multiply_mapping_tRNAfamiles.txt") or die("Cannot open output2 file");

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
			push @{ $StartStop{$splitLine[1]}{$tRNAname} }, $splitLine[8].",".$splitLine[9];
			push @{ $allele{$splitLine[1]}{$tRNAname} }, $splitLine[12];
			push @{ $seq{$splitLine[1]}{$tRNAname} }, $splitLine[13];
        }
    }
}

my $randCounter = 0;
my @duplicatedtRNA;
my %uniqueduptRNA;
my $uniqueCounter = 0;

foreach my $readID (keys %StartStop){
	foreach my $tRNAz (keys %{$StartStop{$readID}}){
		if( $#{$StartStop{$readID}{$tRNAz}} >=1 ){
			push @duplicatedtRNA, $tRNAz;
			$randCounter++;
		}
	}
}

print "Reads with two tRNAs: $randCounter\n";

foreach(@duplicatedtRNA){
	$uniqueduptRNA{$_} = 0;
}

print "Duplicated tRNAs: ";
print scalar keys %uniqueduptRNA;


close(inp0); #Finish reading in BLAST file
