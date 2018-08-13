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
my ($findir,$foutdir,$fshelldir,$fmiddledir,$Nr,$Nt,$Swissprot,$TrEMBL,$Cog,$Kog,$Kegg);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$findir,
				"middir:s"=>\$fmiddledir,
				"od:s"=>\$foutdir,
				"shdir:s"=>\$fshelldir,
				"nr"=>\$Nr,
				"nt"=>\$Nt,
				"swissprot"=>\$Swissprot,
				"trembl"=>\$TrEMBL,
				"cog"=>\$Cog,
				"kog"=>\$Kog,
				"kegg"=>\$Kegg,
				) or &USAGE;
&USAGE unless ($findir and $foutdir and $fmiddledir);

###### paras check ###################
$findir = Cwd::realpath($findir);
$foutdir = Cwd::realpath($foutdir);
$fmiddledir = Cwd::realpath($fmiddledir);

if (defined $fshelldir) {
	$fshelldir = Cwd::realpath($fshelldir);
	mkdir($fshelldir) if (!-d $fshelldir);
}else{
	$fshelldir = getcwd;
}

$findir =~ s/\/$//;
my $Cpu = 50;
my $notename = ` hostname `;
chomp $notename;

my $dirname = dirname($findir);

######get the the best hit ############
mkdir("$dirname/work_sh") if (!-d "$dirname/work_sh");
mkdir($foutdir) if (!-d $foutdir);
mkdir($fmiddledir) if (!-d $fmiddledir);
my $blast_parser = "$Bin/blast_parser.pl";
open TAB1,">$fshelldir/Covertm0_2Tabbest.sh" || die "$!";
open TAB2,">$fshelldir/Covertm0_2Tab.sh" || die $!;

my @blastm0files=glob("$findir/*.blast");
my $seq_file_name = basename("$blastm0files[0]");

$seq_file_name =~ s/\.(\d+?)\.fa\.\S+//;
foreach my $blastfile (sort @blastm0files) {
	my $basename = basename($blastfile);
	if($blastfile=~/NT.blast/){
		print TAB1 "perl $Bin/blast_parser_blastn.pl -nohead -tophit 1 -m 0 -topmatch 1 $blastfile > $fmiddledir/$basename.tab.best && \n";
		print TAB2 "perl $Bin/blast_parser_blastn.pl -nohead -tophit 50 -m 0 -topmatch 1 $blastfile > $fmiddledir/$basename.tab &&\n" ;
	}else{
		print TAB1 "perl $blast_parser -nohead -tophit 1 -m 0 -topmatch 1 $blastfile > $fmiddledir/$basename.tab.best &&\n" ;
		print TAB2 "perl $blast_parser -nohead -tophit 50 -m 0 -topmatch 1 $blastfile > $fmiddledir/$basename.tab &&\n" ;
	}
}
close TAB1;
close TAB2;

&Cut_shell_qsub("$fshelldir/Covertm0_2Tabbest.sh",$Cpu,"1G","scr.q,all.q");
&Cut_shell_qsub("$fshelldir/Covertm0_2Tab.sh",$Cpu,"1G","scr.q,all.q");

############## Convert format tab #########
my @blastm7files=glob("$findir/*.blast.xml");
$seq_file_name = basename("$blastm7files[0]");
$seq_file_name =~ s/\.(\d+?)\.fa\.\S+//;

##get the the best hit
open TAB3,">$fshelldir/Covertm7_2Tabbest.sh" || die "$!";
open TAB4,">$fshelldir/Covertm7_2Tab.sh" || die $!;

foreach my $blastfile (sort @blastm7files) {
	my $outprefix=basename($blastfile);
	$outprefix=~s/\.xml$//;
	print TAB3 "perl $blast_parser -nohead -tophit 1 -topmatch 1 -m 7 $blastfile > $fmiddledir/$outprefix.tab.best &&\n" ;
	print TAB4 "perl $blast_parser -nohead -tophit 20 -topmatch 1 -m 7 $blastfile > $fmiddledir/$outprefix.tab &&\n" ;
}
close TAB3;
close TAB4;

&Cut_shell_qsub("$fshelldir/Covertm7_2Tabbest.sh",$Cpu,"1G","scr.q,all.q");
&Cut_shell_qsub("$fshelldir/Covertm7_2Tab.sh",$Cpu,"1G","scr.q,all.q");

############ Cat Tab OUT #####################
my @files=();
if (defined $Nr) {
	system("cat $fmiddledir/*.NR.blast.tab.best >$foutdir/$seq_file_name.NR.blast.tab.best");
	system("cat $fmiddledir/*.NR.blast.tab >$foutdir/$seq_file_name.NR.blast.tab");
}
if(defined $Nt){
	system("cat $fmiddledir/*.NT.blast.tab.best >$foutdir/$seq_file_name.NT.blast.tab.best");
	system("cat $fmiddledir/*.NT.blast.tab >$foutdir/$seq_file_name.NT.blast.tab");
}
if(defined $Swissprot){
	@files=glob("$fmiddledir/*.Swissprot.blast.tab.best");
	if (scalar @files) {
		system("cat $fmiddledir/*.Swissprot.blast.tab.best >$foutdir/$seq_file_name.Swissprot.blast.tab.best");
		system("cat $fmiddledir/*.Swissprot.blast.tab >$foutdir/$seq_file_name.Swissprot.blast.tab");
	}else{
		system("cat $findir/*.Swissprot.blast.tab >$foutdir/$seq_file_name.Swissprot.blast.tab");
		system("perl -ane 'print if ++\$hash{ \$F[0] }<2' $foutdir/$seq_file_name.Swissprot.blast.tab >$foutdir/$seq_file_name.Swissprot.blast.tab.best");
	}
}
if(defined $TrEMBL ){
	@files=glob("$fmiddledir/*.TrEMBL.blast.tab.best");
	if (scalar @files) {
		system("cat $fmiddledir/*.TrEMBL.blast.tab.best >$foutdir/$seq_file_name.TrEMBL.blast.tab.best");
		system("cat $fmiddledir/*.TrEMBL.blast.tab >$foutdir/$seq_file_name.TrEMBL.blast.tab");
	}else{
		system("cat $findir/*.TrEMBL.blast.tab >$foutdir/$seq_file_name.TrEMBL.blast.tab");
		system("perl -ane 'print if ++\$hash{ \$F[0] }<2' $foutdir/$seq_file_name.TrEMBL.blast.tab >$foutdir/$seq_file_name.TrEMBL.blast.tab.best");
	}
}
if(defined $Cog ){
	@files=glob("$fmiddledir/*.COG.blast.tab.best");
	if (scalar @files) {
		system("cat $fmiddledir/*.COG.blast.tab.best >$foutdir/$seq_file_name.COG.blast.tab.best");
		system("cat $fmiddledir/*.COG.blast.tab >$foutdir/$seq_file_name.COG.blast.tab");
	}else{
		system("cat $findir/*.COG.blast.tab >$foutdir/$seq_file_name.COG.blast.tab");
		system("perl -ane 'print if ++\$hash{ \$F[0] }<2' $foutdir/$seq_file_name.COG.blast.tab >$foutdir/$seq_file_name.COG.blast.tab.best");
	}
}
if(defined $Kog ){
	@files=glob("$fmiddledir/*.KOG.blast.tab.best");
	if (scalar @files) {
		system("cat $fmiddledir/*.KOG.blast.tab.best >$foutdir/$seq_file_name.KOG.blast.tab.best");
		system("cat $fmiddledir/*.KOG.blast.tab >$foutdir/$seq_file_name.KOG.blast.tab");
	}else{
		system("cat $findir/*.KOG.blast.tab >$foutdir/$seq_file_name.KOG.blast.tab");
		system("perl -ane 'print if ++\$hash{ \$F[0] }<2' $foutdir/$seq_file_name.KOG.blast.tab >$foutdir/$seq_file_name.KOG.blast.tab.best");
	}
}
if(defined $Kegg ){
#	system("cat $fmiddledir/*.KEGG.blast.tab.best >$foutdir/$seq_file_name.KEGG.blast.tab.best");
#	system("cat $fmiddledir/*.KEGG.blast.tab >$foutdir/$seq_file_name.KEGG.blast.tab");
	system("cat $findir/*.KEGG.blast.tab >$foutdir/$seq_file_name.KEGG.blast.tab");
}




#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = system "wc -l $shell";
	if ($line<=1000) {
		if ($notename=~/login-0-0/) {
			system "$Bin/qsub-sge.pl $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
		else
		{
			system "ssh login-0-0 $Bin/qsub-sge.pl $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
	if ($line>=1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div = "";
		my $div_index=1;
		my $line_num=0;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=0;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			if ($notename=~/login-0-0/) {
				system "$Bin/qsub-sge.pl $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "ssh login-0-0 $Bin/qsub-sge.pl $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
		}
	}
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
		-id       <dir>  input blast result dir                              required
		-od       <dir>  output dir                                          required
		-middir   <dir>  dir where put the middle result file            required
		-shdir    <dir>  dir where put shell script  default[./]        optional
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
