# Writer:         chaijingchao
# # Program Date:   2016.06.27
#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw ($Bin);
my $usage =
"
Usage:

  Options:

  -h          Help
  -in         classify_summary.txt file
  -out        the output prefix

";

my($in, $out, $help, $info, $Num_filter, $Num_nfl, $Num_fl, $Num_fl_nc, $Num_fl_c);
GetOptions(
          "in=s" => \$in,
          "out=s" => \$out,
          "help|h=s" => \$help
);
if(!$in || $help){
        die "$usage\n";
}
$out ||= "./classify";

open IN, "<$in";
while(<IN>){
             chomp;
             if(/Number of filtered short reads=(\d+)/){
                       $Num_filter = $1;
             }elsif(/Number of non-full-length reads=(\d+)/){
                       $Num_nfl = $1;
             }elsif(/Number of full-length reads=(\d+)/){
                       $Num_fl = $1;
             }elsif(/Number of full-length non-chimeric reads=(\d+)/){
                       $Num_fl_nc = $1;
             }
}       
$Num_fl_c = $Num_fl - $Num_fl_nc;
close IN;
system("/share/public/software/R-3.3.3/bin/Rscript $Bin/classify_pie_plot $Num_filter $Num_nfl $Num_fl_nc $Num_fl_c $out");

