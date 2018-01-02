#!/usr/bin/perl -w
use strict;
my $file=$ARGV[0];
my @array;
#my $line;
my $count=0;
open IN,"$file";
open OUT1 ,">$file.txt";
while(<IN>){
#$line=$_;
if($_=~/^>/){$count++;@array=split /\s/,$_;print OUT1 ">".$array[7]." ".$array[3]." ".$array[4]." ".$array[5]."\n";
}
else {print OUT1 $_;
}
}
close IN;
close OUT1;
