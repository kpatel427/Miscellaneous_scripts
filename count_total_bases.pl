#!/usr/bin/env perl  
#Khushbu Patel 04-27-18
#This script provides number of bases in a multiline fasta file. 

use warnings;
use strict;
use Bio::SeqIO;

#get command-line arguments, or die with a usage statement
my $usage         = "Usage: ./count_total_bases.pl infile |\n infile= filename\n";
my $infile        = shift or die $usage;  #type filename
 
#create one SeqIO object to read in
my $seq_in = Bio::SeqIO->new(
                             -file   => "<$infile",
                             -format => 'fasta',
                             );

							 
	
	my $id;        #scalar to hold sequence identifier    
	my $string;    #scalar to hold the sequence
	my $length =0;
	   	
#write each entry in the input file to STDOUT, while the condition is true
#seq_in is the object that holds the file
#next seq passes each contig to the sequence object $seq
while (my $seq = $seq_in->next_seq ()) 
{ 
	
	#Initialize variables; start the counter at 0 for each contig that comes through the loop  
	
	$id= $seq->id; #Get the header Id of each contig
	$length +=$seq->length;  #Declare variable for finding the length of each contig sequence
	
	$string = $seq->seq;   #Each contig seq is treated as a string
	
	
}


print "Number of total bases = $length\n";