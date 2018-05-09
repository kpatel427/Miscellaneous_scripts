#!/usr/bin/env perl  
# Khushbu Patel 3-16-18
# hqSNP diagnose script / QC lyve-SET/1.1.4f
# Run this script from your main Lyve-SET folder

use warnings;
use strict;
use Bio::SeqIO;



#----------------------------------------------------------------------------------------------------------------------------------------------------#

my $usage = "\nUsage: ./hqSNP_diagnose.pl ./reference/reference.fasta\n\nPlease set the working directory to the main Lyve-SET folder and supply a reference genome!\n\nMake sure you have loaded modules - samtools/1.4.1 and perl/5.16.1-MT!\n";
my $reference = $ARGV[0] or die $usage;

               
# variable declarations

               my $basename;
               my $mapped_reads;
               my $total_mapped_reads;
               my $mapped_reads_lt_30;
               my $not_mapped;
               my $tot_reads;
               my $perc_mapped_reads_gt_30;
               my $perc_mapped_reads;
               my $file;
               my $line_count;
               my $proper_paired;
               my $perc_proper_paired;
               my $total_ref_bases;
			   my $str5;
			   my $perc_reference_coverage;
			   my $average_coverage;
			   my $vcf_out;
               


#.....Setting the paths to fetch files!.........#
my $path_cleaned_reads = "../reads";
my $vcf_path = "../vcf";
my $out_aln_path = "../msa/out.aln.fasta";
my $out_info_path = "../msa/out.informative.fasta";





#Get total number of bases in the Reference

sub gsize()  #Subroutine to calculate size of Genome
{
               #create one SeqIO object to read in
               my $seq_in = Bio::SeqIO->new(
                             -file   => "<$reference",
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

               $total_ref_bases = $length;
               print "\nGenome Size: $total_ref_bases\n";


               }

               
#Printing the reference size         
gsize();


chdir("./bam/");

 #-------------------------------------------------------------------------------------------------------------------------------#
               sub masking;      #function declaration to prevent error: "too early to prototype"
               
               
               
               print "\n\n\n================================================== MASKING INFORMATION ==============================================================\n";
               
               print "\n\n#......Overall masking.....#\n";
               print "-------------------------------------\n";
               
               masking($out_aln_path);              #calling subroutine masking and passing out.aln file as a parameter
               
               print "\n\n#......Masking in informative sites.....#\n";
               print "-------------------------------------\n";

               masking($out_info_path);             #calling subroutine masking and passing out.info file as a parameter


			   
			   
#-------------------------------------------------------------------------------------------------------------------------------#			   
			   

my @files = glob("*.bam");

foreach my $i(@files)
{

               
               ($basename = $i) =~ s/.fastq*.+//;            #Get the basename of files
			   
		       print ("\n------------------------------------------------------------------------------------------------------------------------------------------------------\n");
               print "\n\n\n---> Processing file: $basename...\n";
               print ("\n\n=============================== QC on BAM File =========================\n");
               
               $total_mapped_reads = `samtools view -c -F 4 $i`; #Fetching total number of mapped reads  #NOTE: saving system() output to perl variable, returns 0; it sets the exit status to perl variable; Hence use ``
               $mapped_reads = `samtools view -q 30 $i | wc -l`;  #Fetching number of mapped reads with quality higher than 30
               $mapped_reads_lt_30 = $total_mapped_reads - $mapped_reads; #Number of reads < Q30
               
               
               
               #Calculating total number of reads from clean_reads files
               my $file = "$path_cleaned_reads"."/"."$basename*gz";          
               my $line_count= `zcat $file | wc -l`;   #gettting total number of lines in fastq files
               $tot_reads = $line_count/4;          #Dividing total number of lines by 4, to get total number of reads
               $not_mapped = $tot_reads - $total_mapped_reads; # Calculating number of reads that are not mapped
               $perc_mapped_reads = ($total_mapped_reads/$total_mapped_reads) * 100; #Calculating percentage of mapped reads
               $perc_mapped_reads_gt_30 = ($mapped_reads/$total_mapped_reads) * 100; #Calculating percentage of mapped reads > Q30
			   
			   
			 
               print "\nTotal number of reads in File -> $basename: $tot_reads\n";  # Printing total number of reads
               print "Total number of mapped reads: $total_mapped_reads";
               print "Number of Reads Mapped > Q30: $mapped_reads";
               print "Number of Reads mapped < Q30: $mapped_reads_lt_30\n";
               print "Number of Reads not mapped: $not_mapped\n";
               print "Percent mapped reads: $perc_mapped_reads\n";
               print "Percent mapped reads > Q30: $perc_mapped_reads_gt_30\n";
              

               my $str1 = `samtools stats  $i | head -37 | grep 'reads properly paired:'`;
               $str1 =~s/^SN\s+//;        #Formatting output
               $str1 =~s/#.*//;
               print $str1;                         #reads properly paired:  ###### is stored as one string. Inorder to get the number, Splitting at ':'
               
               my @temp = split ":", $str1; 
               $proper_paired = $temp[1];        #second element of the array is the number
               $proper_paired =~s/\s+//;           #Removing the white spaces before the number
               $perc_proper_paired = ($proper_paired/$mapped_reads) * 100;
               

               
               my $str2 = `samtools stats  $i | head -37 | grep 'inward oriented pairs:'`;
               $str2 =~s/^SN\s+//;                        #Formatting output
               print $str2;
               
               my $str3 = `samtools stats  $i | head -37 | grep 'outward oriented pairs:'`;
               $str3 =~s/^SN\s+//;                        #Formatting output
               print $str3;
               
               my $str4 = `samtools stats  $i | head -37 | grep 'pairs with other orientation:'`;
               $str4 =~s/^SN\s+//;                        #Formatting output
               print $str4;
               
               
               print "Percentage of properly paired reads = $perc_proper_paired\n";
              
               
               $str5 = `samtools stats $i | head -37 | grep 'insert size average:'`;
               $str5 =~s/^SN\s+//;                        #Formatting output
               print $str5;
               
               
               my $baseslt20 = `samtools depth $i| awk '( (\$3 >= 0) && (\$3 <= 19) )'| wc -l`;
               print "Number of bases with mapping coverage less than 20: $baseslt20";

               
               $perc_reference_coverage = ((($total_ref_bases - $baseslt20)/$total_ref_bases)*100);
               print "Percent reference coverage = $perc_reference_coverage\n";
               
      
               
               my $tot_coverage_bases = `samtools depth $i | cut -f 3| awk '{sum += \$_} END {print sum}'`;  # Total coverage of all bases
               print "Total coverage of bases (based off 20x at each position) = $tot_coverage_bases";
               
               $average_coverage = $tot_coverage_bases/$total_ref_bases;
               print "Average Coverage = $average_coverage";


               print "\n\n\n=============================== QC on VCF File =========================\n";
               my $vcf_filename= "$vcf_path"."/".$basename."*.gz";
               
               
               print("Filters applied for sites that do not pass -
a. DP20 = Depth is less than 20, the user-set coverage threshold
b. RF0.95 = Reference variant consensus is less than 0.95, the user-set threshold
c. AF 0.95 = Allele variant consensus is less than 0.95, the user-set threshold
d. isIndel = Indels are not used for analysis in Lyve-SET
e. masked = This site was masked using a bed file or other means
f. str10 = Less than 10% or more than 90% of variant supporting reads on one strand
g. indelError = Likey artifcat due to indel reads at this position\n\n");
               
               
				$vcf_out =`zcat $vcf_filename | grep -v '^#' | cut -f 7 | sort | uniq -c`;
				print("$vcf_out");
				

               }


			  
               
               
               
               
               
               #subroutine to calculate masking from out.aln and out.informative file
               sub masking()
               {
                              my ($file) = @_;
                              
               #create one SeqIO object to read in
               my $seq_in = Bio::SeqIO->new(
                             -file   => "<$file",
                             -format => "fasta",
                             );

                                                                                                         
               my @unsorted1;  #array to hold scalar variables below 
               my $id;        #scalar to hold sequence identifier    
               my $string;    #scalar to hold the sequence
               my $full_id;   #scalar to hold the full id including length etc. 
                              
               #write each entry in the input file to STDOUT, while the condition is true
               #seq_in is the object that holds the file
               #next seq passes each contig to the sequence object $seq
               while (my $seq = $seq_in->next_seq ()) 
               { 
               
                              #Initialize variables; start the counter at 0 for each contig that comes through the loop
               
                              my $count=0;  #total count
                              my $countN=0; #number of ambiguous bases (N)
               
                              $id= $seq->id; #Get the header Id of each contig
                              my $length=$seq->length;  #Declare variable for finding the length of each contig sequence
               
                              $string = $seq->seq;   #Each contig seq is treated as a string
                              $countN = $string =~ s/(N)//sgi;
                              my $percentN = ($countN / $length) *100;
               
                              $full_id= "$id %Ambiguous_bases(N)= $percentN Ambiguous_base_count(N)= $countN Contig_length= $length";
               
                              push (@unsorted1, [$full_id, $length]); #Push each scalar into the unsorted array (identifiers, sequence and percent ambiguous bases) 
                                 #array elements 0,1,2    #Note $num  must also be added separately from the full_id for sorting on GC
                                                                                                                                                                        
               
               }
    #Declare a new array for the sorting the unsorted array. Use $$ to reference the other array. 
               #Sort the array within an array on $num (the second element in the unsorted array) in ascending order (a <=> b) on GC content 
               #my @sorted= sort {$$a[2] <=> $$b[2]} @unsorted;
               
               #Loop through the array.. Print the contig sequences in order of GC content (lowest to highest)

               for (my $i=0;$i < @unsorted1; $i++) 
               {
               
                              my $seqstruct=$unsorted1[$i];   #The current iteration through the array
                              my $full_id=$$seqstruct[0];  #full id is the 0th element in the array 
                              my $sequence=$$seqstruct[1]; #sequence is the 1st element in the array 
               
                              print ">$full_id\n"; 

                              
               }              
               
               
               
               
} #end of masking subroutine
