#!/usr/bin/env perl  
# Khushbu Patel 3-16-18
# This script provides number of ambiguous bases, ratio and percent ambiguous bases for each read in a fastq file.


use warnings;
use strict;

#get command-line arguments, or die with a usage statement
my $usage         = "Usage: perl test1.pl infile \n";
my $infile        = shift or die $usage;  #takes input file


#to create a fasta file from a fastq file
#my $cmd = "cat PNUSAE012394-WAPHL-M4796-180226.fastq.gz.cleaned.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\\t' '\n' > test.fasta";
#system("$cmd");

#print "Fasta file created!\n";

open(FH,$infile) or die "Cannot open file!:$!\n";

my @result;


while(<FH>)
{

	#variables initialization
	my $ID;
	my $countN = 0;
	my $count_tot = 0;
	my $full_ID;
	my $ratio_N = 0;
	my $percent_N = 0;
	my $count_notN = 0;
	
	chomp $_; #removing trailing new lines
	
	
	if($_ =~/^>/) #Lines begining with > to retrieve reads ID 
	{
		$ID = $_; #Fetching Reads ID
		
		push (@result,$ID);  #pushing ID into result array
		
		
		}
		
	else
	{
	
		my @count = $_ =~/n/gi;
		$countN = scalar@count;  #Counting the number of Ns
		push (@result,"Number of Ambiguous bases [Ns]:$countN");  #pushing Number of N's in result array
		
		
		my @chars = split //, $_;
		$count_tot = scalar@chars;   #Counting the total number of bases
		push (@result,"Total number of bases:$count_tot");   #pushing the total number of bases into result array
		
	
		$ratio_N = (($countN)/($count_tot-$countN));   #calculating the ratio of N over everything but N.
		$percent_N = ($countN/$count_tot )*100;   #Calculating percentage of Ns
		
		
		
		push (@result, "Ratio of Ambiguous bases:$ratio_N","Percentage of Ambiguous bases:$percent_N","\n");  #pushing ratio_N and percent_N calculated in the above steps, into result array
	
		}
		
		
		
	
	}
	

	
	#printing the result array
	
	foreach (@result)
	{
		print "$_\n";
		}

	