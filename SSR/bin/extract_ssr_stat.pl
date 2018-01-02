#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
die "Usage: perl $0 <ssr.statistics> <prefix>\n" unless @ARGV==2;

my $statistics=shift;
my $prefix=shift;

my $dir=dirname($statistics);

my $f1=0;my $f2=0; my $f3=0;my $f4=0;
my $basic_stat="$dir/$prefix.ssr_summary.txt";
my $category="$dir/$prefix.ssr_classification.txt";
my $ssr_freq1="$dir/$prefix.ssr_motif_frequency.txt";
my $ssr_freq2="$dir/$prefix.repeat_type_frequency.txt";
open BASIC,">$basic_stat";
open SUMMARY,">$category";
open PATTERN1,">$ssr_freq1";
open PATTERN2,">$ssr_freq2";
my $size; my $n=0;
open STAT,$statistics;
while(<STAT>){
	chomp;
	if(/^RESULTS OF MICROSATELLITE SEARCH/){
		$f1=1;
		$f2=$f3=$f4=0;
		next;
	}
	elsif(/^Distribution to different repeat type classes/){
		$f2=1;
		$f1=$f3=$f4=0;
		next;
	}
	elsif(/^Frequency of identified SSR motifs/){
		$f3=1;
		$f1=$f2=$f4=0;
		next;
	}
	elsif(/^Frequency of classified repeat types/){
		$f4=1;
		$f1=$f2=$f3=0;
		next;
	}elsif(/^Total size of examined sequences \(bp\)\:\s+(\d+)/){
		$size=$1;
	}
	if($size and $n==0){
		print SUMMARY "## Total size (bp):\t$size\n";
		$n=1;
	}
	unless(/^[1-9|a-zA-Z]/){
		next;
	}
	
	if($f1){
		my @array=split /:\s+/;
		print BASIC $array[0]."\t".$array[1]."\n";
	}elsif($f2){
		my @array=split /\t/;
		print SUMMARY $array[0]."\t".$array[1]."\n";
	}elsif($f3){
		print PATTERN1 $_."\n";
	}elsif($f4){
		print PATTERN2 $_."\n";
	}
}
close STAT;
close BASIC;
close SUMMARY;
close PATTERN1;
close PATTERN2;

open(SSR,"$ssr_freq1");
my (@start,@end,@start_pos,@end_pos);
my $head=<SSR>;
chomp $head;
my @head=split /\t/,$head;
my $arrar_num=4; #four type as a group
for (my $i=1;$i+$arrar_num-1<$#head;$i+=$arrar_num){
#	print $i+1,"\t",$head[$i],"\t",$i+$arrar_num,"\t",$head[$i+$arrar_num-1],"\n";
	push @start,$head[$i];
	push @start_pos,$i;
	push @end,$head[$i+$arrar_num-1];
	push @end_pos,$i+$arrar_num-1;
}
$end[$#end]=$head[$#head-1];
$end_pos[$#end_pos]=$#head-1;

my @repeat_name;
my @lengend_name;
for my $i(0..$#start){
	my $tmp="repeat(".$start[$i]."-".$end[$i].")";
	my $leng=$start[$i]."-".$end[$i];
	push @repeat_name,$tmp;
	push @lengend_name,$leng;
}
my $lengend_name=join ",",@lengend_name;
$lengend_name="\"".$lengend_name."\"";
#print $lengend_name,"\n";

my %hash;
while(<SSR>){
	chomp;
	my @tmp=split /\t/, $_;
	my $num=length($tmp[0]);
	for my $i(0..$#repeat_name){
		for my $j($start_pos[$i]..$end_pos[$i]){
			if(defined $tmp[$j] && $tmp[$j] ne "-" && $tmp[$j] ne ""){
				$hash{$num}{$repeat_name[$i]}+=$tmp[$j];
			}
			else{
				$hash{$num}{$repeat_name[$i]}+=0;
			}
		}
	}
}
close SSR;

my @arrra=("Mono-","Di-","Tri-","Tetra-","Penta-","Hexa-","Hepta-","Octo-","Nona-","Deca-");
open(OUT,">$dir/SSR_summary_plot.txt");
print OUT "number\t",join("\t",@repeat_name),"\n";
foreach my $num(sort keys %hash){
	my @tmp;
	for my $i(0..$#repeat_name){
		push @tmp,$hash{$num}{$repeat_name[$i]};
	}
	print OUT $arrra[$num-1],"\t",join("\t",@tmp),"\n";
}

close OUT;









