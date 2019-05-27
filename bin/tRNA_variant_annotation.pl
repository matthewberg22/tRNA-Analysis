#!/usr/bin/perl 
use strict;
use warnings;

#By Patrick O'Donoghue and Matt Berg
#Department of Biochemistry, University of Western Ontario
#May 2019

#Script annotates genomic coordinates with canonical tRNA numbering
#To run enter tRNA_variant_annotation.pl GtRNAdb.txt tRNAstructurenumber.txt nonredundant_mutants.txt
#GtRNAdb.txt contains fasta file of all tRNAs, downloaded from http://gtrnadb.ucsc.edu/GtRNAdb2/index.html in September 2017
#tRNAstructurenumber comes from individual_tRNA_numberingscript.pl
#nonredundant_mutants.txt comes from Allele_Frequencies.pl

#Reads in input files
my $GtRNAdb = $ARGV[0]; 
my $numbering = $ARGV[1];
my $variants = $ARGV[2];

#OUTPUT FILES
system("rm Annotated_nonredundant_mutants.txt"); 
open(out1, ">Annotated_nonredundant_mutants.txt") or die("Cannot open output1 file");

#Reads in GtRNAdb and save tRNA information;
open(inp0, "$GtRNAdb") or die("Cannot open database files: $GtRNAdb, usage is: tRNAvariantScript.pl GtRNAdb.txt tRNAnumbering.txt variants.txt");
   my $genei=0;
   my %chr;
   my %tRNAUpCoord;
   my %tRNADnCoord;
   my %tRNAStrand;
   
while(<inp0>){
	chomp $_;
	  
	my $carrotTest = substr($_,0,1); #Tests if the line starts with > signifying an info line
      
	if ($carrotTest eq ">") {
		my @splitLine = split(" ", $_);
		my @temptRNA = split("_", $splitLine[0]);
		
		my $tRNA = $temptRNA[2];

		#Record Chromosome number for each tRNA gene
        my @tempChr = split(":", $splitLine[4]);
        $chr{$tRNA} = $tempChr[0];

		#Record coordinates for each tRNA gene
        my @tempCoord = split("-", $tempChr[1]);
        $tRNADnCoord{$tRNA} = $tempCoord[0];
        $tRNAUpCoord{$tRNA} = $tempCoord[1];
		
		#Record positive or negative strand
		$tRNAStrand{$tRNA} = $splitLine[5];
         
		
	}
}

close(inp0);



#Reads in numbering file and saves how each tRNA is numbered 5' -> 3';
open(inp1, "$numbering") or die("Cannot open database file $numbering, usage is: tRNAvariantScript.pl GtRNAdb.txt tRNAnumbering.txt variants.txt");
my %numbering;
my $numberingcounter = 0;

while(<inp1>){
    my @splitLine = split " "; #Splits on each space for a line and saves in the array
	my $tRNA = $splitLine[0];
	shift @splitLine;
	
	$numbering{$tRNA} = \@splitLine;
}

close(inp1);


my %tRNAnumbering; #hash variable to store information tRNA information for each scaffold position
my @sequenceline; #array variable to hold the sequence of each tRNA, with each position in the array being one nucleotide
my $seqcounter; #counts through the positions of the tRNA, when information is being assigned to each scaffold position
my $position; #variable to store the scaffold position
my $i; #counter to keep track of what position in the tRNA is being assigned information and when the end of the tRNA has been reached
my @tRNApositionNumber;

foreach my $keys (keys %chr){
	
	foreach my $numberings (keys %numbering){

		if($keys eq $numberings){

				#if the tRNA gene is on the (-) strand, reverses sequence 
				#because even though tRNA sequences are reporter in the 5' to 3' direction of the tRNA, 
				#the scaffold number will be reversed for genes on the - strand
				if($tRNAStrand{$keys} eq '(-)'){
					my @array = @{$numbering{$keys}};
					@tRNApositionNumber = reverse @array;
				}
				
				else{
					@tRNApositionNumber = @{$numbering{$keys}};
				}

			$seqcounter = 0; #sequence counter counts through the numbering array for the tRNA (set to 0 at the start of each tRNA that is being worked)
			$position = $tRNADnCoord{$keys}; #the first position of the tRNA, when this value is reached, every position in that one tRNA has been assigned information
			
			#loops through all the positions in the tRNA
			for ($i = $position; $position <= $tRNAUpCoord{$keys}; $position++){
				my $test = "$chr{$keys}:$position"; #stores the scaffold position to a variable
				$tRNAnumbering{$test} = $tRNApositionNumber[$seqcounter]; #assigns the tRNA number to the position in the scaffold
				$seqcounter++;
			}
		}
	}
}


#Read in list of tRNA gene variants
open(inp2, "$variants") or die("Cannot open database files: $variants, usage is: tRNAvariantScript.pl GtRNAdb.txt tRNAnumbering.txt variants.txt");
	my $vari = 0; #keeps track of which variant is being recorded
	my %variantChr; #hash to extract the chromosome the variant is on
	my %variantCoord; #hash to extract the position of the variant
	my %mut; #hash to store the mutation
	my %refBase; #hash to store the reference bases
	my %mutBase; #hash to store the mutant base
	
print out1 "#tRNA\tMutation\tAF\tSeq\tNumber_Mutations\ttRNA_Position\ttRNA_location\n";

while(<inp2>){
   my $line = $_; #stores the whole line in the variable
   chomp($line); #removes any unwanted/hidden characters at the end of the line
   
   my $check = substr($line, 0, 1);
   
   if($check ne "H"){
		my @splitLine = split("\t", $line);
		
		my @multiplemutations = split(" ", $splitLine[1]);
		
		if(scalar(@multiplemutations) == 1){
	   
		    $vari++; #advances the counter noting the number of variants
		 
		    #Record chromosome of each variant
			my @varChr = split(":", $splitLine[1]);
			my $variant = "$varChr[0]:$varChr[1]";
			my $location;
			
			$location = "No information";
			$location = "01_Acceptor Stem" if "1 2 3 4 5 6 7 66 67 68 69 70 71 72" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "02_D-stem" if "10 11 12 13 22 23 24 25" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "03_D-loop" if "14 15 16 17 17A 17B 17C 18 19 20 20A 20B 20C 21" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "04_Anticodon Stem" if "27 28 29 30 31 39 40 41 42 43" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "05_Anticodon Loop" if "32 33 37 38" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "06_Anticodon" if "34 35 36" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "07_Variable Arm" if "45 46 47" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "07_Variable Arm" if $tRNAnumbering{$variant} =~ m/e/;
			$location = "08_T-Stem" if "49 50 51 52 53 61 62 63 64 65" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "09_T-Loop" if "54 55 56 57 58 59 60" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "10_Discriminator Base" if "73" =~ /\b$tRNAnumbering{$variant}\b/;
			$location = "Intron" if $tRNAnumbering{$variant} =~ m/i/;
		    
			chomp $_;
		   print out1 "$_\t1\t$tRNAnumbering{$variant}\t$location\n";
		   
		}
		
		elsif(scalar(@multiplemutations) > 1){
			my $multiplevariants = "";
			my $multiplelocations = "";
			my $howmanymutants = scalar @multiplemutations;
			
			foreach(@multiplemutations){
				$vari++; #advances the counter noting the number of variants
			 
				#Record chromosome of each variant
				my @varChr = split(":", $_);
				my $variant = "$varChr[0]:$varChr[1]";
				
				my $location;
			
				$location = "No information";
				$location = "01_Acceptor Stem" if "1 2 3 4 5 6 7 66 67 68 69 70 71 72" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "02_D-stem" if "10 11 12 13 22 23 24 25" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "03_D-loop" if "14 15 16 17 17A 17B 17C 18 19 20 20A 20B 20C 21" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "04_Anticodon Stem" if "27 28 29 30 31 39 40 41 42 43" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "05_Anticodon Loop" if "32 33 37 38" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "06_Anticodon" if "34 35 36" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "07_Variable Arm" if "45 46 47" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "07_Variable Arm" if $tRNAnumbering{$variant} =~ m/e/;
				$location = "08_T-Stem" if "49 50 51 52 53 61 62 63 64 65" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "09_T-Loop" if "54 55 56 57 58 59 60" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "10_Discriminator Base" if "73" =~ /\b$tRNAnumbering{$variant}\b/;
				$location = "Intron" if $tRNAnumbering{$variant} =~ m/i/;
				
				$multiplevariants = $multiplevariants."$tRNAnumbering{$variant} ";
				$multiplelocations = $multiplelocations."$location ";
			
			}
			chomp $_;
			print out1 "$_\t$howmanymutants\t$multiplevariants\t$multiplelocations\n";
		}

	}
}

close(inp2);

   
