#!/usr/bin/perl -w
use strict;
use FindBin '$Bin';


if (@ARGV<3){
	die "$0 Trinity.fasta ssr.results outputfile\n";
}
my ($unigene,$ssr_result,$out) = @ARGV;
my ($length,%ssr) = ();

open GENE,"$unigene" || die "$!";
while(<GENE>){
	chomp;
	next if (/^>/);
	my $tmp_length = length $_;
	$length += $tmp_length;
}
close GENE;

open SSR,"$ssr_result" || die "$!";
<SSR>;
while(<SSR>){
	chomp;
	my @line = split /\s+/;
	$line[2] =~ s/\*//;
	$ssr{$line[2]} += 1;
}
close SSR;

open OUT,">$out" || die "$!";
foreach my $type(sort keys %ssr){
	my $ratio = (1000000*$ssr{$type})/$length;
	$ratio = sprintf("%.2f",$ratio);
	print OUT "$type\t$ratio\n";
}
close OUT;

system("Rscript $Bin/draw_ssr_density.R $out");

