#!/usr/bin/perl -w
use strict;
use Getopt::Long;

`export PATH=/share/public/software/bowtie2/bowtie2-2.1.0:\$PATH`;
print "export PATH=/share/public/software/bowtie2/bowtie2-2.1.0:\$PATH\n";

my $dir_assembly_file;
my $dir_reads_file;
my $name_sample = "sample";  #default:sample
my $dir_output ;
my $name_reference = "reference";  #default:reference
my $number_of_processors = 8;  #default:8
my $phred_quality_standard = 33;  #default:33
my $read_type = "p";  #default:pair-end
my $tmp;

my $help;
my $ss=0.5;

GetOptions("h|help" => \$help,
           "a=s" => \$dir_assembly_file,
           "t=s" => \$read_type,
           "r=s" => \$dir_reads_file,
           "o=s" => \$dir_output,
           "n=s" => \$name_reference,
           "s=s" => \$name_sample,
           "p=i" => \$number_of_processors,
             "ss=f"=> \$ss,
	   "tmp=s"=>\$tmp,
           "q=i" => \$phred_quality_standard) or die "Cannot get the input: $!";
	     

if (!defined($dir_assembly_file) || !defined($dir_reads_file) || defined($help))
{
	&usage;
	exit 0;
}

$phred_quality_standard = "--phred$phred_quality_standard-quals";

$dir_output ||= "./$name_sample.RSEM.out";
$dir_output = "./".$dir_output unless ($dir_output =~ /\//);
$dir_output =~ s/\/$//;
mkdir $dir_output unless -d $dir_output;
#print $dir_output;

if($read_type eq "p")
{
  $dir_reads_file =~ s/,/ /;
}

my $dir_bowtie2 = "/share/public/software/bowtie2/bowtie2-2.1.0";
my $dir_rsem_calculate_expression = "/share/public/software/trinity/trinityrnaseq_r20140717/trinity-plugins/rsem/rsem-calculate-expression";
my $dir_temp = $tmp;
if( !-d $dir_temp){die "tmp dir is not exist";}
#my $dir_output_of_extract_transcript_to_gene_map_from_trinity = "$dir_temp/$name_reference.togene";
my $dir_output_of_rsem_prepare_reference = "$dir_temp/$name_reference";
if($read_type eq "s"){
	print "perl $dir_rsem_calculate_expression --bowtie2 --bowtie2-path $dir_bowtie2 $phred_quality_standard -p $number_of_processors --forward-prob $ss --time  $dir_reads_file $dir_output_of_rsem_prepare_reference $dir_output/$name_sample --no-bam-output\n"; 
	!system "perl $dir_rsem_calculate_expression --bowtie2 --bowtie2-path $dir_bowtie2 $phred_quality_standard -p $number_of_processors --forward-prob $ss --time  $dir_reads_file $dir_output_of_rsem_prepare_reference $dir_output/$name_sample --no-bam-output" or die "Something went wrong!";
}elsif($read_type eq "p" ){
	print "perl $dir_rsem_calculate_expression --bowtie2 --bowtie2-path $dir_bowtie2 $phred_quality_standard -p $number_of_processors --forward-prob $ss  --time --paired-end $dir_reads_file $dir_output_of_rsem_prepare_reference $dir_output/$name_sample --no-bam-output\n"; 
	!system "perl $dir_rsem_calculate_expression --bowtie2 --bowtie2-path $dir_bowtie2 $phred_quality_standard -p $number_of_processors --forward-prob $ss --time --paired-end $dir_reads_file $dir_output_of_rsem_prepare_reference $dir_output/$name_sample --no-bam-output" or die "Something went wrong!";
}else{
	print "No type matched! You are expected to enter p for paired-end, s for single-end."
}

sub usage
{
	print STDERR <<USAGE;
===============================================================================
QUANTIFICATION OF TRANSCRIPTS ABUNDANCE BY RSEM
Version: 1.2.15-modified by Fangbinbin 2014-8-6
Use bowtie2 instead of bowtie 
Function: To get readcounts and RPKM for each transcript
Author: Chen Xi
--------------------------------------------------------
USAGE: runRSEM_revised.pl [options]
Options:
	-h/--help		help information

    -a	<string>	Trinity.fasta
	
    -t	<string>	type of reads files, s for single-end 
                    and p for pair-end; default: p
    
    -r	<string>	reads file, *.fastq(if single-end), 
                    *.L.fastq,*.R.fastq(if pair-end)
	
    -o	<string>	output directory, 
                    default:"./Sample_name.RSEM.out";
	
    -n	<string>	name of reference
	
    -s	<string>	name of sample
	
    -q	<int>		phred quality standard, choices: 33 or 64
                    default: 33
	
    -p	<int>		number of processors, default:2 

     -ss <f>		Probability of generating a read from the forward strand of a transcript.
			Set to 1.0 for a strand-specific protocol where all (upstream) reads are
			derived from the forward strand, 0.0 for a strand-specific protocol where 
			all (upstream) read are derived from the reverse strand, 
			or 0.5 for a non-strand-specific protocol. (Default: 0.5)
    -tmp <string>	temp dir of rsem
------------------------------------------------------------------------------------------------------------------
USAGE
}
#----------------------------------------------------------


###############################################################################
################ Graph the Density Picture ####################################
###############################################################################

print "$dir_output\n";

#to prepare readcount and RPKM for RPKM density Graph##

#open IN1, "< $dir_assembly_file";
open IN, "< $dir_output/$name_sample.isoforms.results";
open OUT, "> $dir_output/$name_sample.Readcount_FPKM.xls";
print OUT "transcript_id\tSample_name\tRead_count\tFPKM\n";
<IN>;
while(<IN>){
	chomp;
	my @temp2 = split /\t/;
	print OUT "$temp2[0]\t$name_sample\t$temp2[4]\t$temp2[6]\n"; 
}
close IN;
close OUT;
#R function for RPKM density Graph#########################################
my $densityGraphR =<< "END";
#=======================================================================
setwd("$dir_output")
RPKM_density <- function(dat,logMode = TRUE){
    dat <- read.table(dat, header = TRUE)
    if(logMode) dat\$FPKM<-dat\$FPKM
    .libPaths("/share/public/software/R-3.3.3/lib64/R/library")
    library("ggplot2")
    p <- ggplot(dat)
    if(logMode) {
        p<- p+geom_density(aes(x=log10(FPKM),group=Sample_name,color=Sample_name,fill=Sample_name),alpha=I(1/3))
    }else{
        p<- p+geom_density(aes(x=FPKM,group=Sample_name,color=Sample_name,fill=Sample_name),alpha=I(1/3))
    }
    p<- p + labs(y="Density", title="FPKM density distribution")
    p<- p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
	p <- p + theme(
		panel.background = element_rect(fill = "transparent",colour =NA),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		plot.background = element_rect(fill  = "transparent",colour =NA),
		axis.line=element_line()
	)
    p
}
END

#------------------------------------------------------------------------
open OUT, "> $dir_output/graph.R";
print OUT "$densityGraphR\n";
print OUT "pdf(\"$name_sample.FPKM_density_distribution.pdf\")\n";
print OUT "RPKM_density(dat=\"$name_sample.Readcount_FPKM.xls\")\n";
print OUT "dev.off()\n";
print OUT "png(\"$name_sample.FPKM_density_distribution.png\", type=\"cairo-png\")\n";
print OUT "RPKM_density(dat=\"$name_sample.Readcount_FPKM.xls\")\n";
print OUT "dev.off()\n";
close OUT;
!system "/share/public/software/R-3.3.3/bin/R CMD BATCH $dir_output/graph.R" or die "Something went wrong when graphing the density picture!";
#`samtools view $dir_output/$name_sample.transcript.bam>$dir_output/$name_sample.transcript.sam`;
#`samtools depth $dir_output/$name_sample.transcript.sorted.bam>$dir_output/$name_sample.depth`;


