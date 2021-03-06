#!/usr/bin/perl

=head1 Name

  cog_parser.pl  --  extract the cog number.

=head1 Description

  This program is designed for extracting the cog number after blast.

=head1 Version

  Author: sunjuan, sunjuan@genomics.org.cn
  Version: 1.3,  Date: 2008-1-31

=head1 Usage
	
  perl cog_parser.pl <blast_tab>
  --nohead      do not show the first instruction line.
  --verbose     output verbose information to screen  
  --help        output help information to screen  

=head1 Exmple

  perl cog_parser.pl PLASMID.fasta.ori.glimmer3.pep.blast.tab

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
my $whog = "/share/software/database/COG/whog";
my $fun = "/share/software/database/COG/fun.txt";
my $seq_name = $1 if ($blast_tab=~/([\w\.]+)\.COG\.blast/);
my $outdir = dirname($blast_tab);
#print "$seq_name\n";

my (%COG,%fun);
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
	my ($class,$cog_num,$cog_anno) = ($1,$2,$3) if (/(\[\S+\])\s+(COG\d+)\s+(.+?)\n/s);
	$cog_anno =~ s/\s+/ /g;
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
			$COG{$protein[$j]} = "$org\t$cog_num\t$cog_anno\t$class\t$fun{$class}{category}\t$fun{$class}{defination}";
		}
	}
}
$/="\n";
close WHOG;

##read the tab file and create a file including cog info
open IN,$blast_tab || die "fail $blast_tab";
open OUT,">$outdir/$seq_name.COG.class" || die "fail $seq_name.COG.class";
print OUT "#GeneID\tCOG_protein_NAME\tE_value\tIdentity\tScore\tOrganism\tCOG_id\tCOG_class_defination\tFunction_code\tFunctional_categories\tFunction_class_defination\n" unless (defined $Nohead);
while (<IN>) {
	my @t = split /\t/;
	print OUT "$t[0]\t$t[4]\t$t[13]\t$t[8]\t$t[12]\t$COG{$t[4]}\n" if (exists $COG{$t[4]});
}
close OUT;
close IN;

