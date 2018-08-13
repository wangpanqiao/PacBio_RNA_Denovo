#!/usr/bin/perl

=head1 Name

  kog_parser.pl  --  extract the kog number.

=head1 Description

  This program is designed for extracting the kog number after blast.

=head1 Version

  Author: sunjuan, sunjuan@genomics.org.cn
  Version: 1.3,  Date: 2008-1-31

=head1 Usage
	
  perl kog_parser.pl <blast_tab>
  --nohead      do not show the first instruction line.
  --verbose     output verbose information to screen  
  --help        output help information to screen  

=head1 Exmple

  perl kog_parser.pl PLASMID.fasta.ori.glimmer3.pep.blast.tab

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Nohead,$Verbose,$Help);
GetOptions(
	"nohead"=>\$Nohead,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $blast_tab = shift;
my $whog = "/share/work2/liuying/condatabase/KOG/wkog";
my $fun = "/share/work2/liuying/condatabase/KOG/fun.txt";
my $seq_name = $1 if ($blast_tab=~/([\w\.]+)\.KOG\.blast/);
my $outdir = dirname($blast_tab);
#print "$seq_name\n";

my (%KOG,%fun);
my $category;

##read fun.txt
open FUN,$fun || die "fail $fun";
while (<FUN>) {
	if (/(\[\w\])\s(.+?)\n$/) {
		my ($class,$function) = ($1,$2);
		$fun{$class}{defination} = $function;
		$fun{$class}{category} = $category;
	}else {
		chomp $_;
		$category = $_;
	}
}
close FUN;

##read whog
open WHOG,$whog || die "fail $whog";
$/="_______";
while (<WHOG>) {
	my ($class,$kog_num,$kog_anno) = ($1,$2,$3) if (/(\[\S+\])\s+(KOG\d+)\s+(.+?)\n/s);
	$kog_anno =~ s/\s+/ /g;
	$_ =~ s/.+?\[.+?\n//s;
	my @line = split (/\n/,$_);
	shift @line;
	my $org;
	for (my $i=0;$line[$i];$i++) {
		$org = $1 if ($line[$i]=~/(\S+)\:/);
		$line[$i] =~ s/.+?\://;
		$line[$i] =~ s/^\s+//;
		my @protein = split (/ /,$line[$i]);
		for (my $j=0;$protein[$j];$j++) {
			$fun{$class}{category} = "--" unless (exists $fun{$class}{category});
			$fun{$class}{defination} = "--" unless (exists $fun{$class}{defination});
			$KOG{$protein[$j]} = "$org\t$kog_num\t$kog_anno\t$class\t$fun{$class}{category}\t$fun{$class}{defination}";
		}
	}
}
$/="\n";
close WHOG;

##read the tab file and create a file including kog info
open IN,$blast_tab || die "fail $blast_tab";
open OUT,">$outdir/$seq_name.KOG.class" || die "fail $seq_name.Kog.class";
print OUT "#GeneID\tKOG_protein_NAME\tE_value\tIdentity\tScore\tOrganism\tKOG_ID\tKOG_class_defination\tFunction_code\tFunctional_categories\tFunction_class_defination\n" unless (defined $Nohead);
while (<IN>) {
	chomp;
	s/\s+$//;
	next if (/^$/);
	my @t = split /\t/,$_;
	#print OUT "$t[0]\t$t[4]\t$t[13]\t$t[8]\t$t[12]\t$KOG{$t[1]}\n" if (exists $KOG{$t[1]}); #$KOG{$t[4]} => $KOG{$t[1]}
	print OUT "$t[0]\t$t[1]\t$t[10]\t$t[2]\t$t[11]\t$KOG{$t[1]}\n" if (exists $KOG{$t[1]});
}
close OUT;
close IN;

