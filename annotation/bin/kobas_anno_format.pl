#!/usr/bin/perl -w
use strict;
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
my ($anno,$foutdir,$fkey);
GetOptions(
				"help|h|?" =>\&USAGE,
				"ann:s"=>\$anno,
				"od:s"=>\$foutdir,
				"key:s"=>\$fkey,
				) or &USAGE;
&USAGE unless ($anno and $foutdir and $fkey);
$anno = Cwd::realpath($anno);
$foutdir = Cwd::realpath($foutdir);
mkdir $foutdir if (!-d $foutdir);

my %globemap = (
		"ko01100" => 1,
		"ko01110" => 1,
		"ko01120" => 1,
		);

########database file############
#my $keggsmap = "/share/nas2/database/kegg/pathway/ko";
my $pathclass = "/share/public/database/KEGG/pathway_menu.htext";
my %path_class;
&load_pathclass_ ($pathclass,\%path_class);

######## format out put xls ####
my $gene2ko_file = "$foutdir/$fkey.KEGG.gene2ko.xls";
my $gene2path_file = "$foutdir/$fkey.KEGG.gene2path.xls";

my %path_stat;
my %pathway_stat;
open IN,"$anno" || die $!;
my $part = 0;
$/="--------------------";
while (<IN>) {
	chomp;
	s/^\s+//;s/\s+$//;
	next if (/^$/);
	$part ++;

	if ($part == 1) {
		&sub_gene2ko_ ($_,$gene2ko_file);
	}
	if ($part == 2) {
		&sub_gene2path_ ($_,$gene2path_file,\%path_class,\%path_stat,\%pathway_stat);
	}
	if ($part > 2) {
		print "kobas anno file is wrong! please Check $anno...\n";
		die;
	}
}
close IN;

######### stat path class numbers #########
open OUT,">$foutdir/$fkey.KEGG.pathstat.xls" || die $!;
foreach my $level1 (sort keys %path_stat) {
	print OUT ">$level1 (level1)\n";
	foreach my $level2 (sort keys %{$path_stat{$level1}}) {
		my $genes_num = keys %{$path_stat{$level1}{$level2}};
		print OUT "\t$level2 (level2)\t$genes_num\n";
		foreach my $kos (sort keys %{$pathway_stat{$level1}{$level2}}) {
			print OUT "\t\t$kos\t$path_class{$kos}{name}\t$pathway_stat{$level1}{$level2}{$kos}\n";
		}
	}
}
close OUT;

################pick up kegg map##############



#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub load_pathclass_ {
	my ($file,$class) = @_;

	open IN,"$file" || die $!;
	my ($level1,$level2);
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;
		next if (/^$/ || /^\#/);
		if (/^A /) {
			$level1 = (split /\s+/,$_,2)[1];
		}
		if (/^B /) {
			$level2 = (split /\s+/,$_,2)[1];
		}
		if (/^C /) {
			my ($id,$def) = (split /\s+/,$_,3)[1,2];
			$id = "ko"."$id";
			$class -> {$id}{level1} = $level1;
			$class -> {$id}{level2} = $level2;
			$class -> {$id}{name} = $def;
		}
	}
}

sub sub_gene2ko_ {
	my ($info,$gene2ko_out) = @_;

	my %gene2ko;
	my @kos = split /\n+/,$info;
	foreach my $info_line (@kos) {
		$info_line =~ s/^\s+//;
		$info_line =~ s/\s+$//;
		next if ($info_line =~ /^$/ || $info_line =~ /^\#/);
		my ($gene,$ko_info) = split /\t+/,$info_line;
		next if ($ko_info =~ /None$/);
		$ko_info =~ s/\|http:/\thttp:/;
		$gene2ko{$gene} .= "$ko_info\t";
	}

	open OUT,">$gene2ko_out" || die $!;
	print OUT "#GeneID\tKEGG_entry\tKEGG_url\n";
	foreach my $g (sort keys %gene2ko) {
		$gene2ko{$g} =~ s/\s+$//;
		print OUT "$g\t$gene2ko{$g}\n";
	}
	close OUT;
}

sub sub_gene2path_ {
	my ($info,$gene2path_out,$path_clas,$stat,$pathstat) = @_;

	my %path;
	my @infos = split/\/\/\/\//;
	foreach (@infos) {
		s/^\s+//;s/\s+$//;
		next if (/^$/);
		my @tmp = split /\n+/,$_;
		next if ($#tmp <= 1);
		my ($gene,$ko);
		$gene = (split /\t+/,$tmp[0],2)[1];
		$ko = (split /\t+/,$tmp[1],2)[1];
		$gene =~ s/^\s+//; $gene =~ s/\s.*$//;
		$ko =~ s/^\s+//; $ko =~ s/\s.*$//;

		for (2..$#tmp) {
			if ($tmp[$_] =~ /^Pathway:/) {
				$tmp[$_] =~ s/^Pathway://;
			}
			$tmp[$_] =~ s/^\s+//;
			next if ($tmp[$_] !~ /KEGG PATHWAY/);
			my ($path_name,$path_acc) = (split /\t+/,$tmp[$_])[0,2];
			next if (exists $globemap{$path_acc});
			if (defined $path_clas -> {$path_acc}) {
				$stat -> {$path_clas -> {$path_acc}{level1}}{$path_clas -> {$path_acc}{level2}}{$gene} = 1;
				$pathstat -> {$path_clas -> {$path_acc}{level1}}{$path_clas -> {$path_acc}{level2}}{$path_acc} ++;
			}
			$path{$path_acc}{gene} .= "$gene+";
			$path{$path_acc}{K} .= "$ko+";
			$path{$path_acc}{name} = $path_name;
			$path{$path_acc}{number} ++;
		}
	}

	open OUT,">$gene2path_out" || die $!;
	print OUT "#Pathway\tpathway_ko\tGenes_Number\tGeneIDs\tKEGG_entrys\n";
	foreach my $path_acc (sort keys %path) {
		$path{$path_acc}{gene} =~ s/\+$//;
		$path{$path_acc}{K} =~ s/\+$//;
		print OUT "$path{$path_acc}{name}\t$path_acc\t$path{$path_acc}{number}\t$path{$path_acc}{gene}\t$path{$path_acc}{K}\n";
	}
	close OUT;
}

sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[2014/12/31]
	Contact:fred routine <fred_routine\@163.com.cn>
	Description:	Extract KEGG Pathway and KOs file From kobas anno file;

	Usage:
		-ann         The blast_tab Out of seq with KEGG           must be given

		-od          Output dir                                   must be given

		-key         Prefix of OUT files (eg. key.path key.ko)    must be given

		-h          Help document
USAGE
	print $usage;
	exit;
}
