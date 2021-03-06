#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2012 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2012 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2012 
my $version="1.0.0";
my $BEGIN=time();

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info
# ------------------------------------------------------------------
my ($Database_cfg,$odir,$LOG,$VERBOSE,$species,$HELP);
my ($nr,$nt,$Swissprot,$TrEMBL,$GO,$COG,$KOG,$Kegg,$all);
my ($path,$Blastp,$Blast_cpu,$Blast_e,$Blast_cut);
my ($sh_dir,$Result_dir,$Tab_dir,$div_dir,$step);

my @anno;

GetOptions(
		"Database_cfg:s"=>\$Database_cfg,
		"cog"=>\$COG,
		"kog"=>\$KOG,
		"all"=>\$all,
		"nr"=>\$nr,
		"nt"=>\$nt,
		"swissprot"=>\$Swissprot,
		"trembl"=>\$TrEMBL,
		"kegg"=>\$Kegg,
		"GO"=>\$GO,
		"path=s"=>\$path,
		"Blastp:s"=>\$Blastp,
		"Blast_cpu:i"=>\$Blast_cpu,
		"Blast_e:f"=>\$Blast_e,
		"Blast_cut:i"=>\$Blast_cut,
		"step:s"=>\$step,
		"od:s"=>\$odir,
		"log:s"=>\$LOG,
		"verbose"=>\$VERBOSE,
                "species:s"=>\$species,
		"help"=>\$HELP
	) or &USAGE;

&USAGE if (!$path || !$odir || !$species || $HELP) ;

$step||=1;

###------------------软件路径----------------------------###
# all the program are in this dir                    
chomp $Bin;

#=================== 一些参数 ================================
&MKDIR($odir);
$odir=&ABSOLUTE_DIR($odir);
$Database_cfg=&ABSOLUTE_DIR($Database_cfg);


$Blast_cpu = $Blast_cpu // 20;
$Blastp = $Blastp // "blastx";
$Blast_e = $Blast_e // 1e-5;
$Blast_cut = $Blast_cut // 200;


#my @arr = <$path/Transcripts.fa>;
#die "$path/Transcripts.fa not exists...\n" if (@arr < 1);
#my $Query = $arr[0];
#my $Q_name=basename $Query;
my $Query = "$path/UniIso.fa";
unless (-e $Query){die "$path/UniIso.fa not exists...\n";}
my $Q_name=basename $Query;

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n$programe_dir Start Time :[$Time_Start]\n\n";

#######日志文件
if (defined $LOG) {
	open (LOG, ">>$LOG") ||  die "Can't write $LOG: $!\n";
}
else{
	open (LOG, ">>LOG.txt") ||  die "Can't write LOG: $!\n";
}


&MKDIR("$odir/Result");
$Result_dir="$odir/Result";
&MKDIR("$odir/02.gene-annotation");
$Tab_dir="$odir/02.gene-annotation";
&MKDIR("$odir/Query_Seq_div");
$div_dir="$odir/Query_Seq_div";
&MKDIR("$div_dir/work_sh");
$sh_dir="$div_dir/work_sh";

################################ Program
my %FA1;
&LOAD_SEQ($Query,\%FA1);

if ($step==1) {
	print "Cut Query fa file...\n";
	&CUTFA(\%FA1,$div_dir,$Blast_cut,$Q_name);
	$step=2;
}

print STDERR "Gene annotation ... \n\n" if defined $VERBOSE;
my $program="02.gene-annotation";

if (( !defined $nr)&& ( !defined $nt) && ( !defined $COG) && (!defined $KOG) && ( !defined $Swissprot) && ( !defined $TrEMBL) && ( !defined $Kegg) && (! defined $GO) && (!defined $all) ) {
	LOGFILE(1,$program);
	print  BOLD RED "$program error: You must choose database to blast\n";
	exit;
}

my $database_choose;
$database_choose .= " --cog"       if (defined $COG);
$database_choose .= " --kog"       if (defined $KOG);
$database_choose .= " --nr"        if (defined $nr || defined $all);
$database_choose .= " --nt"        if (defined $nt || defined $all);
$database_choose .= " --swissprot" if (defined $Swissprot || defined $all);
$database_choose .= " --trembl"    if (defined $TrEMBL || defined $all);
$database_choose .= " --kegg"      if (defined $Kegg || defined $all);

if ($step==2) {
	my $worksh;
	$worksh .= "perl $Bin/bin/blast_database.pl -Q_dir $div_dir -name $Q_name -evalue $Blast_e -p $Blastp -cfg $Database_cfg ";
	$worksh .= "$database_choose";
	$worksh.=" --cpu $Blast_cpu -TAB $Tab_dir\n";

	runcmd ($program,"$worksh");
	$step=3;
}

if ($step==3) {
	runcmd ($program,"perl $Bin/bin/blast_format.pl -id $div_dir/$Q_name.div -od $Tab_dir -shdir $sh_dir  -middir $div_dir/mid $database_choose ");
	$step=4;
}

if ((defined $GO || defined $all) && $step==4) {
	$Q_name =~ /(\w+?)[_\-.].*/;
	my $go_graph_mark = $1;
	runcmd ($program,"perl $Bin/bin/Blast2GO.pl -id $div_dir/$Q_name.div/ -od $Result_dir -sh_dir $sh_dir -k $Q_name");      
        runcmd ($program,"sed -e 's/, /\t/g' $Result_dir/$Q_name.GO.list.xls>$Result_dir/$Q_name.GO.list_tmp.xls");
        runcmd ($program,"perl $Bin/bin/go_classification.pl $Result_dir/$Q_name.GO.list_tmp.xls $Result_dir");
        runcmd ($program,"/share/public/software/R-3.3.3/bin/Rscript $Bin/bin/plot_go_classification_ggplot2.R $Result_dir/GO_classification_count.txt $Result_dir");
	runcmd ($program,"perl $Bin/bin/go_level2_stat.pl -gene2GO $Result_dir/$Q_name.GO.list.xls -lable $go_graph_mark -od $Result_dir -prefix $Q_name -color 3");
	$step=5;
}

if ((defined $Kegg || defined $all) && $step==5) {
	###########Kegg代谢图 ming'an sun
	if (defined $Kegg || defined $all) {
		runcmd ($program,"export PYTHONPATH=/share/public/software/kobas-3.0//src/:/share/public/software/Python-2.7.13/lib/python2.7/site-packages/:/share/public/software/Python-2.7.13/lib/python2.7/site-packages/:/share/public/software/WDLIB/lib/FALCON_unzip_0619.2017/lib/python2.7/site-packages/:/share/public/software/falcon_unzip/lib/python2.7/site-packages/:/share/public/software/kobas-3.0/src/:/share/public/software/falcon_unzip/lib/python2.7/site-packages/ && python2.7 $Bin/bin/annotate.py -i $Tab_dir/$Q_name.KEGG.blast.tab -t blastout:tab -s ko -o $Tab_dir/$Q_name.KEGG.kobas.annotate -r 20 -n 5");
		runcmd ($program,"perl $Bin/bin/kobas_anno_format.pl -ann $Tab_dir/$Q_name.KEGG.kobas.annotate -od $Tab_dir -key $Q_name");
                runcmd ($program,"perl $Bin/bin/kegg_classification_v3.pl $Tab_dir/$Q_name.KEGG.kobas.annotate $Result_dir $species");
	}
	$step=6;
}

if ($step==6 && (defined $COG || defined $KOG)) {
	###########COG分类
	if (defined $COG) {
		runcmd ($program,"perl $Bin/bin/cog_parser.pl $Tab_dir/$Q_name.COG.blast.tab.best");
		runcmd ($program,"perl $Bin/bin/CogFunClassDrawer.pl -i $Tab_dir/$Q_name.COG.class -db COG -o $Tab_dir/$Q_name.COG.class.svg -png");
	}elsif (defined $KOG) {
		runcmd ($program,"perl $Bin/bin/kog_parser.pl $Tab_dir/$Q_name.KOG.blast.tab.best");
		runcmd ($program,"perl $Bin/bin/CogFunClassDrawer.pl -i $Tab_dir/$Q_name.KOG.class -db KOG -o $Tab_dir/$Q_name.KOG.class.svg -png");
	}
	$step=7;
}


########统计基因汇总基因注释结果（NR,COG,KEGG,Swissprot）
if ($step==7) {
	runcmd ($program,"perl $Bin/bin/anno_result_format.pl -id $Tab_dir -od $Result_dir -query $Q_name $database_choose");
	runcmd ($program,"perl $Bin/bin/anno_integrate.pl -gene $Query -id $Result_dir -od $Result_dir -key $Q_name");
	LOGFILE(0,$program);
}

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\n$programe_dir End Time :[$Time_End]\n\n";
&Runtime($BEGIN);

#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+
sub runcmd
{
	my ($program,$cmd)=@_;
	open (SH, ">>work.sh") ||  die "Can't write work: $!\n";
	print SH "$cmd \n";
	system($cmd) && LOGFILE(1,$program); 
	close SH;
}

sub LOGFILE
{
	my $flog=shift;
	my $program=shift;
	my $Time= sub_format_datetime(localtime(time()));
	print LOG "[$Time +0800]\ttask\t0\tstart\t$program\tStart to analysis......\n";
	if($flog==0){
		print LOG "[$Time +0800]\ttask\t0\tend\t$program\tDone.\n";
	}else{
		print LOG "[$Time +0800]\ttask\t0\terror\t$program\tAt least one $program in this section is in error.\n";
		close LOG;
		exit;
	}
}
close LOG;

sub LOAD_SEQ 
{
	my ($fa,$info) = @_;

	open IN,"$fa" || die $!;
	$/='>';
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/ || /^\#/);
		my ($head,$seq)=split/\n+/,$_,2;
		my $id=(split/\s+/,$head)[0];
		$info->{$id}=$seq;
	}
	$/="\n";
	close IN;
}

sub CUTFA 
{
	my ($fa,$od,$cut,$name) = @_;

	&MKDIR("$od/$name.div");
	my %seq=%$fa;
	my @aa=sort(keys %seq);
	my $index=0;
	LAB: for (my $i=1;;) {
		my $num=0;
		open OUT,">$od/$name.div/$name.$i.fa" || die $!;
		for ($index..$#aa) {
			$index++;
			if ($num<$cut) {
				print OUT ">$aa[$_]\n$seq{$aa[$_]}\n";
				$num++;
			}
			if ($num>=$cut) {
				$num=0;
				$i++;
				close OUT;
				if ($index==$#aa+1) {
					last;
				}
				else {
					next LAB;
				}
			}
		}
		if ($num) {
			close OUT;
		}
		last;
	}
}

sub parse_config
{ # load config file
	my $config_file= shift;
	my $DataBase= shift;
	
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/\s+//;s/\s+$//;s/\r$//;
		next if(/$/ or /\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$DataBase->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

sub USAGE 
{
	print <<"	Usage End.";
	Program:Blast with DataBase Annotate denovo Genes
	Version: $version
	Contact: Meng Fei <fred_routine\@163.com>

	Description:

      -cog             search against Cog database (protokaryon)
      -kog             search against Kog database (eukaryon)

      -all             search against all database
            -nr              search against Nr database
            -nt              search against NT database
            -swissprot       search against SwissProt database
            -trembl          search against TrEMBL database
            -kegg            search against Kegg database
            -GO              run Blast2GO stat GO Annotation ( -nr is required )

      -Database_cfg          Annotation Database DIR Config
      -od                    OUT DIR
			-path                  The path of project for search 'Unigene'

      -step                  Program Start step
                1   Start from begin                Default
                2   Start from Blast With Database
                3   Start from Get fixed format file form blast result
                4   Start from Blast2GO
                5   Start from Result tackle
                6   Start from Anno Result Reorganize & Stat

      -log <file>             set log file
      -verbose                show detailed information
      -species                the kegg species, animal,plant,fungi
      -help

	Example:
	  perl  Gene_Func_Anno_Pipline.pl -all -Database_cfg  DataBase_config.txt

	Usage End.
	exit;
}
