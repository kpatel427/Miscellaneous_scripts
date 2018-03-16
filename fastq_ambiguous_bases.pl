#!/usr/bin/env perl  
# Khushbu Patel 3-16-18
# This script provides number of ambiguous bases, ratio and percent ambiguous bases for each read in a fastq file. It generatess an out_report.txt file. It also prompts if user wants to store the fasta file. 


use warnings;
use strict;

#get command-line arguments, or die with a usage statement
my $usage         = "Usage: ./test1.pl infile.fastq \n";
my $infile        = shift or die $usage;  #takes input file


#to create a fasta file from a fastq file
my $cmd = "cat $infile | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\\t' '\n' > test.fasta";
system("$cmd");

#print "Fasta file created!\n";

open(FH,"test.fasta") or die "Cannot open file!:$!\n";

my @result;

#storing the file name without extension.
my $basename;
($basename = $infile) =~ s/\.[^.]+$//;




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
	
	open(FH1, '>', 'out_report.txt');
	
	foreach (@result)
	{
		print FH1 "$_\n";
		}
	
	close(FH1);
	
	print "Done!\n";
	
	
	#Prompts if user wants to save the fasta file
	print "Do you want to save the fasta file? (Y/N) \n";
	my $a = <STDIN>;
	
	if ($a =~/^y|yes|n|no$/i)
    {
        if ($a =~/^no?$/i) 
            { 
                print "File deleted!\n";
				system("rm test.fasta");  #deletes fasta file if user inputs n/N/NO/no
            } 
        if ($a =~/^y(es)?$/i)
            {
                print "file not deleted\n";
				system("mv test.fasta $basename.fasta"); #renaming test.fasta 
            }
    } else {
        #Validates user input and returns message if invalid option is provided
		system("rm test.fasta");
		print "\n\n=====================================\n";
        print "You Have Entered and Incorrect Option\n";
        print "Valid Options are [Y/N]\n";
        print "=====================================\n\n";
        
    }
	
