#!/usr/bin/perl -w
###############################################
#The script is used to add hyperlinks in the k-
#egg pathway enrichment list.Certainly,it also
#can be appled for other lists
################################################
use strict;
use Spreadsheet::WriteExcel;

unless(@ARGV==2){
	print "\n\tUsage:\n\t\tperl $0 <pathway_enrichment.xls> <out.xls>\n";
	print "\n\tExample:\n\t\tperl $0 /share/work2/staff/jiangdezhi/script/RNA_seq/nr_nt_sw_kegg_kog_ann/pipeline/scramble_case_vs_Tpc1115Scrb_ctrl_Down_pathway_enrichment.xls out.xls\n\n";
	exit;
}
 
my ($infile,$outfile)=@ARGV;
my $row=0;
my $col=0;

# Create a new workbook called simple.xls and add a worksheet
my $workbook  = Spreadsheet::WriteExcel->new($outfile);
my $worksheet = $workbook->add_worksheet();

open IN,$infile or die;
while(<IN>){
     chomp;
     $row++;
     my @arr=split /\t/;
     if($row==1){
	$col=@arr;
	# The general syntax is write($row, $column, $token). Note that row and
	# column are zero indexed
	for(my $i=1;$i<=$col;$i++){
		$worksheet->write($row-1,$i-1,$arr[$i-1]);
	}
     }else{
	for(my $i=1;$i<=$col-1;$i++){
		$worksheet->write($row-1,$i-1,$arr[$i-1]);
	}
		#Write hyperlinks
		$worksheet->write($row-1,$col-1,$arr[$col-1]);
    }
}           
close IN;
