#!/usr/bin/perl -w
#use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0";
my @Times = localtime();
#######################################################################################
my $Time_Start = sub_format_datetime(localtime(time()));
print STDOUT "$Script start at:[$Time_Start]\n";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($anno,$foutdir,$fkey,$db,$spec);
GetOptions(
				"help|h|?" =>\&USAGE,
				"anno:s"=>\$anno,
				"od:s"=>\$foutdir,
				"prefix:s"=>\$fkey,
				"db:s"=>\$db,
				"spe"=>\$spec,
				) or &USAGE;
&USAGE unless ($anno and $foutdir and $fkey and $db);

######## paras input ##############
$foutdir = Cwd::realpath($foutdir);
mkdir $foutdir if (!-d $foutdir);
my $nowdir = Cwd::realpath();

######## database file############
my %eval_intra = (
		"1" => [0,1e-150],
		"2" => [1e-150,1e-100],
		"3" => [1e-100,1e-50],
		"4" => [1e-50,1e-5],
		);

my %iden_intra = (
		"1" => [0,50],
		"2" => [50,70],
		"3" => [70,90],
		"4" => [90,95],
		"5" => [95,100],
		);

my %eval;
my %iden;
my %spe if (defined $spec);
my $all_anno;

open IN,"$anno" || die $!;
while (<IN>) {
	chomp;
	s/\s+$//;s/^\s+//;
	next if (/^$/ || /^\#/ || /^Unigene/);
	my @tmp = split/\t+/,$_;
	$all_anno ++;
#	$tmp[3] =~ /\((.*)\)/;
	my $identity = $tmp[3];
	foreach my $index (sort {$a<=>$b} keys %iden_intra) {
		if ($identity > ${$iden_intra{$index}}[0] && $identity <= ${$iden_intra{$index}}[1]) {
			$iden{$index} ++;
		}
	}
	if (defined $spec) {
		my $species;
		if ($db eq "NR") {
			$tmp[5] =~ /\[(\w+\s\w+)/;
			$species = $1;
		}
		if ($db eq "Swissprot") {
			$tmp[5] =~ /OS=(\w+\s\w+)/;
			$species = $1;
		}
		if (!defined $species) {
			$species = "Undef";
		}
		$spe{$species} ++;
	}
	if ($tmp[2] == 0) {
		$eval{"0"} ++;
		next;
	}
	foreach my $index (sort {$a<=>$b} keys %eval_intra) {
		if ($tmp[2] > ${$eval_intra{$index}}[0] && $tmp[2] <= ${$eval_intra{$index}}[1]) {
			$eval{$index} ++;
		}
	}
}
close IN;

open OUT,">$foutdir/$fkey"."_$db"."_Evalue.stat.xls" || die $!;
print OUT "#evalue_intraval\tGenes_number\tPercent\n";
foreach my $index (sort {$a<=>$b} keys %eval) {
	if ($index == 0) {
		print OUT "0\t$eval{0}\t";
		printf OUT "%.2f",$eval{"0"}/$all_anno*100;
		print OUT "%\n";
	}
	else {
		my $str = "${$eval_intra{$index}}[0]"."~"."${$eval_intra{$index}}[1]";
		print OUT "$str\t$eval{$index}\t";
		printf OUT "%.2f",$eval{$index}/$all_anno*100;
		print OUT "%\n";
	}
}
close OUT;

open OUT,">$foutdir/$fkey"."_$db"."_Identity.stat.xls" || die $!;
print OUT "#Identity_intraval(%)\tGenes_number\tPercent\n";
foreach my $index (sort {$a<=>$b} keys %iden) {
	my $str = "${$iden_intra{$index}}[0]"."~"."${$iden_intra{$index}}[1]"."(%)";
	print OUT "$str\t$iden{$index}\t";
	printf OUT "%.2f",$iden{$index}/$all_anno*100;
	print OUT "%\n";
}
close OUT;

if (defined $spec) {
	my %top_stat;
	my $top_index = 1;
	my $Others = 0;
	foreach my $spec (sort {$spe{$b} <=> $spe{$a}} keys %spe) {
		if ($top_index <= 10) {
			$top_stat{$spec} = $spe{$spec};
		}
		else {
			$Others += $spe{$spec};
		}
		$top_index ++;
	}

	open OUT,">$foutdir/$fkey"."_$db"."_Species.stat.xls" || die $!;
	print OUT "#Species\tGenes_number\tPercent\n";
	foreach my $spec (sort {$top_stat{$b} <=> $top_stat{$a}} keys %top_stat) {
		print OUT "$spec\t$top_stat{$spec}\t";
		printf OUT "%.2f",$top_stat{$spec}/$all_anno*100;
		print OUT "%\n";
	}
	if ($Others > 0) {
		print OUT "Others\t$Others\t";
		printf OUT "%.2f",$Others/$all_anno*100;
		print OUT "%\n";
	}
	close OUT;
}


#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {
	my $usage=<<"USAGE";
	Program: $Script
	Version: $version: [2015/1/22]
	Contact: Mengfei <fred_routine\@163.com.cn>
	Description: Stat Blast Evalue & Identity & Species distribution and graph Pie plot;

	Usage:
		-anno     The formated Out of blast with NR or Swissprot         must be given

		-db       graph lable of Pie plot (NR or Swissprot)              must be given

		-od       Output dir                                             must be given

		-prefix   Prefix of OUT files (eg. prefix_Evalue.stat.xls)       must be given

		--spe     stat and graph species distributiuon                   optical

		-h          Help document
USAGE
	print $usage;
	exit;
}
