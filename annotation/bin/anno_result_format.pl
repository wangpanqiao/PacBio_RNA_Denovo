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
my ($findir,$foutdir,$Q_name,$Nr,$Nt,$Swissprot,$TrEMBL,$Cog,$Kog,$Kegg);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$findir,
				"od:s"=>\$foutdir,
				"query:s"=>\$Q_name,
				"nr"=>\$Nr,
				"nt"=>\$Nt,
				"swissprot"=>\$Swissprot,
				"trembl"=>\$TrEMBL,
				"cog"=>\$Cog,
				"kog"=>\$Kog,
				"kegg"=>\$Kegg,
				) or &USAGE;
&USAGE unless ($findir and $foutdir and $Q_name);

###### paras check ###################
$findir = Cwd::realpath($findir);
$findir =~ s/\/$//;
my $dirname = dirname($findir);
$foutdir = Cwd::realpath($foutdir);

############ copy and format results ######
if (-f "$findir/$Q_name.COG.class" || defined $Cog) 
{
	system "cp $findir/$Q_name.COG.class $foutdir/$Q_name.COG_class.xls";
	system "cp $findir/$Q_name.COG.class.svg $foutdir/$Q_name.COG.class.svg";
	system "cp $findir/$Q_name.COG.class.png $foutdir/$Q_name.COG.class.png";
	system "cp $findir/$Q_name.COG.class.stat.xls $foutdir/$Q_name.Cog_class.stat.xls";
#	system "cp $findir/$Q_name.COG.blast.tab $foutdir/$Q_name.COG.blast.xls";
#	system "cp $findir/$Q_name.COG.blast.tab.best $foutdir/$Q_name.COG.blastTophit.xls";
}

if (-f "$findir/$Q_name.KOG.class" || defined $Kog) 
{
	system "cp $findir/$Q_name.KOG.class $foutdir/$Q_name.KOG_class.xls";
	system "cp $findir/$Q_name.KOG.class.svg $foutdir/$Q_name.KOG.class.svg";
	system "cp $findir/$Q_name.KOG.class.png $foutdir/$Q_name.KOG.class.png";
	system "cp $findir/$Q_name.KOG.class.stat.xls $foutdir/$Q_name.KOG_class.stat.xls";
#	system "cp $findir/$Q_name.KOG.blast.tab $foutdir/$Q_name.KOG.blast.xls";
#	system "cp $findir/$Q_name.KOG.blast.tab.best $foutdir/$Q_name.KOG.blastTophit.xls";
}

if (-f "$findir/$Q_name.KEGG.blast.tab.best" || defined $Kegg) {
	system "cp $findir/$Q_name.KEGG.blast.tab $foutdir/$Q_name.KEGG.blast.xls";
	system "cp $findir/$Q_name.KEGG.gene2ko.xls $foutdir/";
	system "cp $findir/$Q_name.KEGG.gene2path.xls $foutdir/";
	system "cp $findir/$Q_name.KEGG.pathstat.xls $foutdir/";
}

############## NR NT TrEMBL Swissprot #######
if (defined $Swissprot) {
	open (SWISS,"$findir/$Q_name.Swissprot.blast.tab")or die "cant open file $findir/$Q_name.Swissprot.blast.tab";
	open (OUT,">$findir/$Q_name.Swissprot.blast.tab.anno")or die "cant open file $findir/$Q_name.Swissprot.blast.tab.anno";
	print OUT "#GeneID\tSwissprot_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<SWISS>)
	{
		chomp;
		&ann;
	}
	close OUT;
	close SWISS;

	open (SWISS,"$findir/$Q_name.Swissprot.blast.tab.best")or die "cant open file $findir/$Q_name.Swissprot.blast.tab.best";
	open (OUT,">$findir/$Q_name.Swissprot.blast.tab.best.anno")or die "cant open file $findir/$Q_name.Swissprot.blast.tab.best.anno";
	print OUT "#GeneID\tSwissprot_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<SWISS>)
	{
		chomp;
		&ann;
	}
	close OUT;
	close SWISS;
	###### annotation of each library ######
	system "cp $findir/$Q_name.Swissprot.blast.tab.anno $foutdir/$Q_name.Swissprot.blast.xls";
	system "cp $findir/$Q_name.Swissprot.blast.tab.best.anno $foutdir/$Q_name.Swissprot.anno.xls";
}

if (defined $TrEMBL) {
	open (TREMBL,"$findir/$Q_name.TrEMBL.blast.tab")or die "cant open file $findir/$Q_name.TrEMBL.blast.tab";
	open (OUT,">$findir/$Q_name.TrEMBL.blast.tab.anno")or die "cant open file $findir/$Q_name.TrEMBL.blast.tab.anno";
	print OUT "#GeneID\tTrEMBL_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<TREMBL>)
	{
		chomp;
		&ann;
	}
	close OUT;
	close TREMBL;

	open (TREMBL,"$findir/$Q_name.TrEMBL.blast.tab.best")or die "cant open file $findir/$Q_name.TrEMBL.blast.tab.best";
	open (OUT,">$findir/$Q_name.TrEMBL.blast.tab.best.anno")or die "cant open file $findir/$Q_name.TrEMBL.blast.tab.best.anno";
	print OUT "#GeneID\tTrEMBL_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<TREMBL>)
	{
		chomp;
		&ann;
	}
	close OUT;
	close TREMBL;
	###### annotation of each library ######
	system "cp $findir/$Q_name.TrEMBL.blast.tab.anno $foutdir/$Q_name.TrEMBL.blast.xls";
	system "cp $findir/$Q_name.TrEMBL.blast.tab.best.anno $foutdir/$Q_name.TrEMBL.anno.xls";
}

if (defined $Nr) {
	open (NR,"$findir/$Q_name.NR.blast.tab")or die "cant open file $findir/$Q_name.NR.blast.tab";
	open (OUT,">$findir/$Q_name.NR.blast.tab.anno")or die "cant open file $findir/$Q_name.NR.blast.tab.anno";
	print OUT "#GeneID\tNR_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<NR>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close NR;

	open (NR,"$findir/$Q_name.NR.blast.tab.best")or die "cant open file $findir/$Q_name.NR.blast.tab.best";
	open (OUT,">$findir/$Q_name.NR.blast.tab.best.anno")or die "cant open file $findir/$Q_name.NR.blast.tab.best.anno";
	print OUT "#GeneID\tNR_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<NR>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close NR;
	###### annotation of each library ######
	system "cp $findir/$Q_name.NR.blast.tab.anno $foutdir/$Q_name.NR.blast.xls";
	system "cp $findir/$Q_name.NR.blast.tab.best.anno $foutdir/$Q_name.NR.anno.xls";
}

if (defined $Nt) {
	open (NT,"$findir/$Q_name.NT.blast.tab")or die "cant open file $findir/$Q_name.NT.blast.tab";
	open (OUT,">$findir/$Q_name.NT.blast.tab.anno")or die "cant open file $findir/$Q_name.NT.blast.tab.anno";
	print OUT "#GeneID\tNT_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<NT>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close NT;

	open (NT,"$findir/$Q_name.NT.blast.tab.best")or die "cant open file $findir/$Q_name.NT.blast.tab.best";
	open (OUT,">$findir/$Q_name.NT.blast.tab.best.anno")or die "cant open file $findir/$Q_name.NT.blast.tab.best.anno";
	print OUT "#GeneID\tNT_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<NT>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close NT;
	###### annotation of each library ######
	system "cp $findir/$Q_name.NT.blast.tab.anno $foutdir/$Q_name.NT.blast.xls";
	system "cp $findir/$Q_name.NT.blast.tab.best.anno $foutdir/$Q_name.NT.anno.xls";
}

#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub ann
{
	my @anno = split /\t/, $_;
	print OUT "$anno[0]\t$anno[4]\t$anno[13]\t$anno[8]\t$anno[12]\t$anno[15]\n";
}

sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[2015/1/20]
	Contact:fred routine <fred_routine\@163.com>
	Descriptions: Get fixed format file according to blast result file.
	Options:
		-id       <tabs_dir>       input blast result dir                              required
		-od       <result_dir>     output dir                                          required
		-query    <query_fa_name>  blast fas name                                      required
		-nr              search against Nr annotation                   optional
		-nt              search against Nt annotation                   optional
		-swissprot       search against SwissProt annotation            optional
		-trembl          search against TrEMBL annotation               optional
		-cog             search against Cog annotation                  optional
		-kog             search against Kog annotation                  optional
		-kegg            search against Kegg annotation                 optional
		-h               Help

USAGE
	print $usage;
	exit;
}
