#!/usr/bin/perl -w
# Author:Shilai Xing 
# Program name: singleline.pl
# Description: creates a single line sequence record after each sequence name line

open (IN,"<$ARGV[0]") || die ("\nError: Couldn't open fasta file of multiple sequence lines!\n\n");

my $filename = $ARGV[0];
#open (IN,"<$filename") || die ("\nError: Couldn't open source file containing original FASTA sequences !\n\n");
#$filename =~ s/\.Unigenes\.fa//;
open (OUT1,">$filename.single.fa");
open (OUT2,">$filename.name.convert");

my ($line,$count);
while (<IN>)
  {
$count++;
$line=$_;
chomp $line;
if($line=~/^>/ && $.==1){print OUT1 ">SSR_".$count."\n"; $line=~s/iso.*//;print OUT2 $line."\t"."SSR_".$count."\n";} 
elsif($line=~/^>/ && $.>1) {print OUT1 "\n".">SSR_".$count."\n"; $line=~s/iso.*//; print OUT2 $line."\t"."SSR_".$count."\n";}
else { print OUT1 $line;} 
}
close IN;
close OUT1;
close OUT2;

