#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
#######################################################################################

my $notename=`hostname`;chomp $notename;

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($od,$Index,$id,$worksh,$cpu,$step);
GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$od,
				"id:s"=>\$id,
				"sh_dir=s"=>\$worksh,
				"cpu=s"=>\$cpu,
				"s=s"=>\$step,
				"k:s"=>\$Index,
				) or &USAGE;
&USAGE unless ($od and $Index and $id);

$step=$step || 1;
$id=&ABSOLUTE_DIR($id);
$cpu=$cpu || 30;
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
$worksh=$worksh || "work_sh";
&MKDIR($worksh);
$worksh=&ABSOLUTE_DIR($worksh);

&MKDIR("$od/Blast2go");
my $blast2go="/share/public/software/blast2go/b2g4pipe/";
my $obo="/share/public/database/GO/go.obo";

my @files=sort glob("$id/*NR.blast.xml");

open (SH,">$worksh/blast2go.sh")||die "$!";
my $j=1;

foreach my $file (@files) {
	print SH "/share/public/software/jdk1.8.0_11/bin/java -cp $blast2go/*:$blast2go/ext/*: es.blast2go.prog.B2GAnnotPipe -in $file -out $od/Blast2go/$Index.$j -prop $blast2go/b2gPipe.properties -v -annot &&\n";
	$j++;
}

close SH;

if ($step==1) {
	&Shell_qsub("$worksh/blast2go.sh","scr.q,all.q",$cpu,$worksh);
	$step=2;
}

if ($step==2) {
	my @stat=sort glob("$od/Blast2go/*.annot");
	foreach my $file (@stat) {
		if($file=~/$Index\.(\d+)\./){
			`cat $file >> "$od/$Index.annot"`;
		}
	}
	$step=3;
}

############################################################################################################����.annot��.obo
if ($step==3) {
	open (AN,"$od/$Index.annot") or die $!;
	my (%Query,%GO);

	while (<AN>) {
		chomp;
		next if (/^$/);
		my ($query,$go_id)=(split/\t/,$_)[0,1];
		$Query{$query}{$go_id}=1;
		$GO{$go_id}{$query}=1;
	}
	close AN;

	open (OBO,$obo)||die "$!";
	$/="[Term]";
	my ($go_id,$go_name,$go_class,$anno);
	my %GO_Info;
	my %GO_anno;
	while(<OBO>){
		chomp;
		next if(/^$/);
		my @Term_info=split /\n+/,$_;
		foreach (@Term_info) {
			if($_=~/^id: (.+)/){
				$go_id=$1;
			}elsif($_=~/^name: (.+)/){
				$go_name=$1;
			}elsif($_=~/^namespace: (.+)/){
				$go_class = $1;
			}elsif($_=~/^def: \"(.+)\"/){
				$anno=$1;
				$GO_Info{$go_id}{CLASS}=$go_class;
				$GO_Info{$go_id}{NAME}=$go_name;
				$GO_Info{$go_id}{ANNO}=$anno;
				$GO_anno{$go_id}="$go_id: $go_name ($go_class);";
			}
		}
	}
	$/="\n";
	close OBO;

	my %GO_stat;
	open OUT1,">$od/$Index.GO.list.xls"||die"$!";
	open OUT2,">$od/$Index.GO.anno.xls"||die "$!";
	print OUT2 "#GeneID\tGOTerms_num\tGO_Anno\n";
	foreach my $gene (sort {$a cmp $b} keys %Query) {
		print OUT1"$gene\t";
		my @go_list=(keys %{$Query{$gene}});
		print OUT2 "$gene\t";
		my $go_id_str;
		my $go_anno_str;
		my $num = 0;
		foreach my $go_id (sort {$a cmp $b} keys %{$Query{$gene}}) {
			if (exists $GO_anno{$go_id}) {
				$go_id_str .= "$go_id, ";
				$go_anno_str .= "$GO_anno{$go_id}\t";
				$GO_stat{$GO_Info{$go_id}{CLASS}}{$go_id}{$gene}=1;
				$num ++;
			}
		}
		$go_id_str =~ s/,\s+$//;
		$go_anno_str =~ s/\s+$//;
		print OUT1 "$go_id_str\n";
		print OUT2 "$num\t$go_anno_str\n";
	}
	close OUT1;
	close OUT2;

	open OUT3,">$od/$Index.GO_tree.stat.xls"||die "$!";
	print OUT3"#GOTerm\tGO_def\tUnigene_number\tUnigene_ID\n";
	foreach my $go_class (sort keys %GO_stat) {
		print OUT3 "$go_class\n";
		foreach my $go_term (sort keys %{$GO_stat{$go_class}}) {
			my @genes=keys %{$GO_stat{$go_class}{$go_term}};
			my $gene_id_str=join ", ",@genes;
			my $num=@genes;
			print OUT3 "$go_term\t$GO_Info{$go_term}{NAME}\t$num\t$gene_id_str\n";
		}
	}

	close (OUT3) ;
}

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub Shell_qsub
{ # &Shell_qsub($sh,$qeue,$cpu,$shdir);
	my $sh = shift;
	my $qeue = shift;
	my $cpu = shift;
	my $shdir = shift;

	if ($notename=~/login-0-0/)
	{
		chdir $shdir;
		system "$Bin/qsub-sge.pl $sh --queue $qeue --reqsub -maxproc $cpu --independent";
	}
	else
	{
		chdir $shdir;
		system "ssh login-0-0 $Bin/qsub-sge.pl $sh --queue $qeue --reqsub --maxproc $cpu --independent" ;
	}
}

################################################################
sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:        $version
Contact:        fred routine <fred_routine\@163.com>
Program Date:   2015.1.20
Description:    This program is used to tackle xml(m7 blast) to GO Annotation;
Usage:
                 Options:
                     -id <dir>      input files directory, blast result in format -m7               forced;
  
                     -od <dir>      output file directory                                           forced;
  
                     -k <str>       out file prefix name                                            forced;

                     -sh_dir <str>  shell dir                                                       default ./work_sh;

                     -s  <int>       Step Start the program
                         1       from Blast2GO, tackle blast m7 Result          default;
                         2       from cat qsub Blast2GO's Reault file;
                         3       from GO Annotation stat ;

                     -h         Help

USAGE
	print $usage;
	exit;
}
