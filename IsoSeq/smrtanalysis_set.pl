#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $usage =
"
Usage:

  Options:

  -h          Help
  -minfull            <int>   a minimum full pass, default is 0
  -minpredit          <int>   minimum predited accuracy, default is 75
  -minseq             <int>   a minimum read length, default is 300
  -cDNA               <str>   cDNA_size,Estimated cDNA size,the choice is [under1k,between1k2k,between2k3k,above3k]
  -hq                 <str>   hq_quiver_min_accuracy,Minimum allowed quiver accuracy to classify an isoform as 
                              hiqh-quality, default is 0.99

";

my($input,$minfull,$minpredit,$minseq,$cDNA,$hq,$output,$help);
GetOptions(
          "input=s" => \$input,
          "minfull=i" => \$minfull,
          "minpredit=i" => \$minpredit,
          "minseq=i" => \$minseq,
          "cDNA=s" => \$cDNA,
          "hq=i" => \$hq,
          "output=s" => \$output,
          "h|help"=>\$help
);
if(!$input || !$cDNA || !$output || $help){
        die "$usage\n";
}

$minfull ||= 0;
$minpredit ||= 75;
$minseq ||= 300;
$hq ||= 0.99;

open INPUT, "<$input";
open OUTPUT, ">$output";

while(<INPUT>){
           chomp;
           if(/parameter1/){
           print OUTPUT "                <value>$minfull</value>\n"
           }elsif(/parameter2/){
           print OUTPUT "                <value>$minpredit</value>\n"
           }elsif(/parameter3/){
           print OUTPUT "                <value>$minseq</value>\n"
           }elsif(/parameter4/){
           print OUTPUT "                <value>$cDNA</value>\n"
           }elsif(/parameter5/){
           print OUTPUT "                <value>$hq</value>\n"
           }else{
           print OUTPUT "$_\n"
           }
}
close INPUT;
close OUTPUT;
