#!/usr/bin/env perl  
# Khushbu Patel 3-16-18
# BAM diagnose script / QC lyve-SET/1.1.4f

use warnings;
use strict;

#----------------------------------------------------------------------------------------------------------------------------------------------------#

my $usage         = "Usage: perl bam_diagnose.pl reference.fasta \n";
my $reference = $ARGV[0];

#Prompts if user wants to load the samtools module
	print "Do you want to load module samtools/1.4.1? (Y/N) Skip this step if you have already loaded the module! \n";
	my $a = <STDIN>;  #type y/yes/Yes/n/no/No
	
	if ($a =~/^y|yes|n|no$/i)
    {
        if ($a =~/^no?$/i) 
            { 
                print "Skipping module loading step...\n";
            } 
        if ($a =~/^y(es)?$/i)
            {
                print "Loading module...";
				system("module load samtools/1.4.1"); #loading samtools module
				
				
            }
    } else {
        #Validates user input and returns message if invalid option is provided
		print "\n\n=====================================\n";
        print "You Have Entered and Incorrect Option\n";
        print "Valid Options are [Y/N]\n";
        print "=====================================\n\n";
		exit;
        
    }

	
	
# variable declarations

	my $basename;
	my $mapped_reads;
	my $tot_reads;
	my $perc_mapped_reads;
	my $file;
	my $line_count;
	


#.....Setting the paths to fetch files!.........#

my $path_cleaned_reads = "~/EDLB/projects/Surveillance/Shigella/Local_WA_flexneri_march2018/cleaned_reads/180312_clusterA_march2018";

#get total number of bases in reference
my $ref_reads = `grep -v ">" $reference | wc | awk '{print \$3-\$1}'`; #NOTE: saving system() output to perl variable, returns 0; it sets the exit status to perl variable; Hence use ``


#Setting file paths
my $vcf_path = "~/EDLB/projects/Surveillance/Shigella/Local_WA_flexneri_march2018/lyve-SET_v1.1.4f_ext-ref_95-3008_clusterA_3-12-18/vcf";

my $out_aln_path = "/scicomp/home/oix2/EDLB/projects/Surveillance/Shigella/Local_WA_flexneri_march2018/lyve-SET_v1.1.4f_ext-ref_95-3008_clusterA_3-12-18/msa/out.aln.fasta";

my $out_info_path = "/scicomp/home/oix2/EDLB/projects/Surveillance/Shigella/Local_WA_flexneri_march2018/lyve-SET_v1.1.4f_ext-ref_95-3008_clusterA_3-12-18/msa/out.informative.fasta";



my @files = glob("*.bam");

foreach my $i(@files)
{

	
	($basename = $i) =~ s/.fastq*.+//;
	
	print "\n\n\n---> Processing file: $basename...\n";
	print "\n\n=============================== QC on BAM File =========================\n";
	
	
	$mapped_reads = `samtools view -q 30 $i | wc -l`;  #Fetching number of mapped reads with quality higher than 30
	
	
	#Calculating total number of reads from clean_reads files
	my $file = "$path_cleaned_reads"."/"."$basename"."*.gz";	
	my $line_count= `zcat $file | wc -l`;   #gettting total number of lines in fastq files
	$tot_reads = $line_count/4;          #Dividing total number of lines by 4, to get total number of reads
	my $not_mapped = $tot_reads - $mapped_reads; # Calculating number of reads that are not mapped
	$perc_mapped_reads = ($mapped_reads/$tot_reads) * 100; # Calculating percentage of mapped reads
	
	
	print "\n#......Mapped reads....#\n";
	print "-------------------------------------\n";
	print "total number of reads in File -> $basename: $tot_reads\n";  # Printing total number of reads
	print "Number of Mapped Q30: $mapped_reads";
	print "Number of reads not mapped: $not_mapped\n";
	print "Percent mapped reads: $perc_mapped_reads\n";
	
	
	print "\n#......Read Pairing.....#\n";
	print "-------------------------------------\n";

	my $str1 = `samtools stats  $i | head -37 | grep 'reads properly paired:'`;
	$str1 =~s/^SN\s+//; 	#Formatting output
	$str1 =~s/#.*//;
	print $str1;
	
	my $str2 = `samtools stats  $i | head -37 | grep 'inward oriented pairs:'`;
	$str2 =~s/^SN\s+//;		#Formatting output
	print $str2;
	
	my $str3 = `samtools stats  $i | head -37 | grep 'outward oriented pairs:'`;
	$str3 =~s/^SN\s+//;		#Formatting output
	print $str3;
	
	my $str4 = `samtools stats  $i | head -37 | grep 'pairs with other orientation:'`;
	$str4 =~s/^SN\s+//;		#Formatting output
	print $str4;
	
	
	
	print "\n#......Insert Size.....#";
	print "\n-------------------------------------\n";
	
	my $str5 = `samtools stats $i | head -37 | grep 'insert size average:'`;
	$str5 =~s/^SN\s+//;		#Formatting output
	print $str5;
	


	print "\n\n#......Reference Coverage.....#\n";
	print "-------------------------------------\n";
	
	my $baseslt20 = `samtools depth $i| awk '( (\$3 >= 0) && (\$3 <= 19) )'| wc -l`;
	print "Number of bases with mapping coverage less than 20: $baseslt20";


	my $total_ref_bases = `awk '/^>/ {if (seqlen){print seqlen}; seqlen=0;next; } { seqlen = seqlen +length(\$0)}END{print seqlen}' $reference`;   #calculating total number of bases in reference.fasta
	my $perc_reference_coverage = ((($total_ref_bases - $baseslt20)/$total_ref_bases)*100);
	print "Percent reference coverage = $perc_reference_coverage";
	
	

	print "\n\n#......Average Coverage.....#\n";
	print "-------------------------------------\n";
	
	my $tot_coverage_bases = `samtools depth $i | cut -f 3| awk '{sum += \$_} END {print sum}'`;  # Total coverage of all bases
	print "Total coverage of bases = $tot_coverage_bases";
	
	my $average_coverage = $tot_coverage_bases/$total_ref_bases;
	print "Average Coverage = $average_coverage";


	print "\n\n\n=============================== QC on VCF File =========================\n";
	my $vcf_filename= "$vcf_path"."/".$basename."*.gz";
	#print "\nProcessing file: $basename...\n\n";
	
	print "Filters applied for sites that do not pass -
a. DP20 = Depth is less than 20, the user-set coverage threshold
b. RF0.95 = Reference variant consensus is less than 0.95, the user-set threshold
c. AF 0.95 = Allele variant consensus is less than 0.95, the user-set threshold
d. isIndel = Indels are not used for analysis in Lyve-SET
e. masked = This site was masked using a bed file or other means
f. str10 = Less than 10% or more than 90% of variant supporting reads on one strand
g. indelError = Likey artifcat due to indel reads at this position\n\n";
	
	
	system("zcat $vcf_filename | grep -v '^#' | cut -f 7 | sort | uniq -c");

	} 


	
	
	print "\n\n\n=============================== MASKING INFORMATION =========================\n";
	
	print "\n\n#......Overall masking.....#\n";
	print "-------------------------------------\n";
	
	my @array;
	my $str5 = `cat /scicomp/home/oix2/EDLB/projects/Surveillance/Shigella/Local_WA_flexneri_march2018/lyve-SET_v1.1.4f_ext-ref_95-3008_clusterA_3-12-18/log/launch_set.log | grep 'reportMaskedGenomes: '`;
	
	# .............Parsing launch_set.log file........#
	@array = split ":|\n", $str5;
			
	foreach(@array)	
	{
			
			if($_ =~/\sPNU.+/g)
			{
			print $_."\n";
			}
	
		
	}		
	
	print "\n\n#......Masking in informative sites.....#\n";
	print "-------------------------------------\n";
	
	
	open(FH1,$out_info_path) or die "Cannot open file!:$!\n";

	my @result;

	while(<FH1>)
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
		
		push (@result,"\n$ID\t");  #pushing ID into result array
		
		
		}
		
	else
	{
	
		my @count = $_ =~/n/gi;
		$countN = scalar@count;  #Counting the number of Ns
		push (@result,"$countN\t");  #pushing Number of N's in result array
		
		
		my @chars = split //, $_;
		$count_tot = scalar@chars;   #Counting the total number of bases
		push (@result,"$count_tot\t");   #pushing the total number of bases into result array
		
	
		$ratio_N = (($countN)/($count_tot-$countN));   #calculating the ratio of N over everything but N.
		$percent_N = ($countN/$count_tot )* 100;   #Calculating percentage of Ns
		
		
		
		push (@result, "$ratio_N\t","$percent_N\t");  #pushing ratio_N and percent_N calculated in the above steps, into result array
	
		}
		
		
	
	}


	
	
	#............... printing Masking information .....................#

print"\t\t\t\t\t\t\tNumber of Ambiguous bases [Ns]:\t Total number of bases:\t\t Ratio of Ambiguous bases:\t Percentage of Ambiguous bases:\t";	
	
foreach(@result)
{
	print("\t$_\t");

	}
  