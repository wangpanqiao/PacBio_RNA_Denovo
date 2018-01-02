#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $ver="1.0.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);

##############################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作；
my %opts;
GetOptions(\%opts,"gene=s","id=s","od=s","key=s","h");
if (!defined($opts{gene}) || !defined($opts{id}) || !defined($opts{od}) || !defined($opts{key}) || defined($opts{h})) {
	&help;
	exit;
}

#######################

my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n$programe_dir Start Time :[$Time_Start]\n\n";

######################
my $gene=&ABSOLUTE_DIR($opts{gene});
my $anno_dir=&ABSOLUTE_DIR($opts{id});
&MKDIR($opts{od});
my $odir=&ABSOLUTE_DIR($opts{od});
my $cds_key=$opts{key};

my %Gene_length;
open IN,"$gene" || die;
$/='>';
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	my ($title,$seq)=split/\n+/,$_,2;
	$seq=~s/\s+//g;
	my $id=(split/\s+/,$title)[0];
	my $length=length $seq;
	$Gene_length{$id}=$length;
}
$/="\n";
close IN;

my %Anno_Gene;
my %Anno_Stat;
my %Anno_head;
my %Venn_list;

my @Anno_files=glob "$anno_dir/*.anno.xls";

my $ko_anno=glob "$anno_dir/*.gene2ko.xls";
push @Anno_files,$ko_anno if (defined $ko_anno);
my %Anno_Base;

foreach my $file (@Anno_files) {
	my $database_key;
	if ($file=~/$cds_key\.KEGG\.gene2ko\.xls/) {
		$database_key="KEGG";
	}
	elsif ($file=~/$cds_key/ && $file=~/anno/) {
		$file=~m/$cds_key\.(.*).anno.xls/;
		$database_key=$1;
	}
	$Anno_Base{$database_key}=1;
	if ($database_key!~/GO/) {
		open IN,"$file" || die;
		my $head = <IN>; chomp $head;
		my @head_info = split/\t+/,$head;
		$Anno_head{$database_key}{str} = join"\t",@head_info[1..$#head_info];
		$Anno_head{$database_key}{num} = $#head_info-1;
		while (<IN>) {
			chomp;
			next if (/^$/ || /^\#/);
			my @annotate=split/\t+/,$_,2;
			my $id=$annotate[0];
			my $anno=$annotate[-1];
			$Venn_list{$database_key}{$id}=1;
			$Anno_Stat{All}{$id}=1;
			$Anno_Stat{$database_key}{Anno}++;
			$Anno_Gene{$id}{$database_key}=$anno;
                        if ($Gene_length{$id}>=0 && $Gene_length{$id}<1000) {
                                $Anno_Stat{$database_key}{0}++;
                        }
                        if ($Gene_length{$id}>=1000 && $Gene_length{$id}<2000) {
                                $Anno_Stat{$database_key}{1000}++;
                        }
                        if ($Gene_length{$id}>=2000 && $Gene_length{$id}<3000) {
                                $Anno_Stat{$database_key}{2000}++;
                        }
                        if ($Gene_length{$id}>=3000 && $Gene_length{$id}<6000) {
                                $Anno_Stat{$database_key}{3000}++;
                        }
                        if ($Gene_length{$id}>=6000) {
                                $Anno_Stat{$database_key}{6000}++;
                        }
		}
		close IN;
	}
	else {
		open IN,"$file" || die;
		my $head = <IN>; chomp $head;
		my @head_info = split/\t+/,$head;
		$Anno_head{$database_key}{str} = join"\t",@head_info[1..$#head_info];
		$Anno_head{$database_key}{num} = $#head_info-1;
		while (<IN>) {
			chomp;
			next if (/^$/ || /^\#/);
			my ($id,$sum,$anno)=split/\t+/,$_,3;
			$Venn_list{$database_key}{$id}=1;
			$Anno_Stat{All}{$id}=1;
			$Anno_Stat{$database_key}{Anno}++;
			$anno=~s/^\d+\t+//;
			$anno=~s/\t+/ /g;
			$Anno_Gene{$id}{$database_key}="$sum\t"."$anno";
                        if ($Gene_length{$id}>=0 && $Gene_length{$id}<1000) {
                                $Anno_Stat{$database_key}{0}++;
                        }
                        if ($Gene_length{$id}>=1000 && $Gene_length{$id}<2000) {
                                $Anno_Stat{$database_key}{1000}++;
                        }
                        if ($Gene_length{$id}>=2000 && $Gene_length{$id}<3000) {
                                $Anno_Stat{$database_key}{2000}++;
                        }
                        if ($Gene_length{$id}>=3000 && $Gene_length{$id}<6000) {
                                $Anno_Stat{$database_key}{3000}++;
                        }
                        if ($Gene_length{$id}>=6000) {
                                $Anno_Stat{$database_key}{6000}++;
                        }
		}
		close IN;
	}
}


if (-f "$anno_dir/$cds_key.COG_class.xls") {
	my $cog="$anno_dir/$cds_key.COG_class.xls";

	$Anno_Base{COG}=1;

	open COG,"$cog" || die;
	while (<COG>) {
		chomp;
		next if (/^$/ || /^\#/) ;
		my @anno=split/\t+/,$_;
		$Venn_list{COG}{$anno[0]}=1;
		$Anno_Stat{All}{$anno[0]}=1;
		$Anno_Stat{COG}{Anno}++;
		$Anno_Gene{$anno[0]}{COG}="$anno[-3]\t$anno[-1]";
                if ($Gene_length{$anno[0]}>=0 && $Gene_length{$anno[0]}<1000) {
                        $Anno_Stat{COG}{0}++;
                }
                if ($Gene_length{$anno[0]}>=1000 && $Gene_length{$anno[0]}<2000) {
                        $Anno_Stat{COG}{1000}++;
                }
                if ($Gene_length{$anno[0]}>=2000 && $Gene_length{$anno[0]}<3000) {
                        $Anno_Stat{COG}{2000}++;
                }
                if ($Gene_length{$anno[0]}>=3000 && $Gene_length{$anno[0]}<6000) {
                        $Anno_Stat{COG}{3000}++;
                }
                if ($Gene_length{$anno[0]}>=6000) {
                        $Anno_Stat{COG}{6000}++;
                }
	}
	close COG;
}

if (-f "$anno_dir/$cds_key.KOG_class.xls") {
	my $kog="$anno_dir/$cds_key.KOG_class.xls";

	$Anno_Base{KOG}=1;

	open KOG,"$kog" || die;
	while (<KOG>) {
		chomp;
		next if (/^$/ || /^\#/) ;
		my @anno=split/\t+/,$_;
		$Venn_list{KOG}{$anno[0]}=1;
		$Anno_Stat{All}{$anno[0]}=1;
		$Anno_Stat{KOG}{Anno}++;
		$Anno_Gene{$anno[0]}{KOG}="$anno[-3]\t$anno[-1]";
                if ($Gene_length{$anno[0]}>=0 && $Gene_length{$anno[0]}<1000) {
                        $Anno_Stat{KOG}{0}++;
                }
                if ($Gene_length{$anno[0]}>=1000 && $Gene_length{$anno[0]}<2000) {
                        $Anno_Stat{KOG}{1000}++;
                }
                if ($Gene_length{$anno[0]}>=2000 && $Gene_length{$anno[0]}<3000) {
                        $Anno_Stat{KOG}{2000}++;
                }
                if ($Gene_length{$anno[0]}>=3000 && $Gene_length{$anno[0]}<6000) {
                        $Anno_Stat{KOG}{3000}++;
                }
                if ($Gene_length{$anno[0]}>=6000) {
                        $Anno_Stat{KOG}{6000}++;
                }
	}
	close KOG;
}

open OUT,">$odir/Integrated_Function.annotation.xls" || die;
print OUT "#GeneID";
foreach (sort keys %Anno_Base) {
	if (/COG/) {
		print OUT "\tCOG_class\tCOG_class_annotation";
		next;
	}
	if (/KOG/) {
		print OUT "\tKOG_class\tKOG_class_annotation";
		next;
	}
	else {
		print OUT "\t$Anno_head{$_}{str}";
	}
}

print OUT "\n";

my @NA_anno = ("--","--","--","--","--","--","--","--","--","--");
foreach (keys %{$Anno_Stat{All}}) {
	$Anno_Stat{zz}{Anno}++;
	my $id=$_;
	print OUT "$id";
	if ($Gene_length{$id}>=0 && $Gene_length{$id}<1000) {
		$Anno_Stat{zz}{0}++;
	}
        if ($Gene_length{$id}>=1000 && $Gene_length{$id}<2000) {
                $Anno_Stat{zz}{1000}++;
        }
        if ($Gene_length{$id}>=2000 && $Gene_length{$id}<3000) {
                $Anno_Stat{zz}{2000}++;
        }
        if ($Gene_length{$id}>=3000 && $Gene_length{$id}<6000) {
                $Anno_Stat{zz}{3000}++;
        }
        if ($Gene_length{$id}>=6000) {
                $Anno_Stat{zz}{6000}++;
        }
	foreach (sort keys %Anno_Base) {
		if (/COG/ && !defined $Anno_Gene{$id}{COG}) {
			print OUT "\t--\t--";
			next;
		}
		if (/KOG/ && !defined $Anno_Gene{$id}{KOG}) {
			print OUT "\t--\t--";
			next;
		}
		elsif (defined $Anno_Gene{$id}{$_}) {
			printf OUT "\t$Anno_Gene{$id}{$_}";
			next;
		}
		else {
			print OUT "\t";
			printf OUT join "\t",@NA_anno[0..$Anno_head{$_}{num}];
			next;
		}
	}
	print OUT "\n";
}
close OUT;

open STAT,">$odir/Function_Annotation.stat.xls" || die;
print STAT "Anno_Database\tAnnotated_Number\t0<=length<1000\t1000<=length<2000\t2000<=length<3000\t3000<=length<6000\tlength>=6000\n";
foreach (sort keys %Anno_Stat) {
	next if (/All/) ;
	if ($_!~/zz/) {
		print STAT "$_"."_"."Annotation\t$Anno_Stat{$_}{Anno}\t$Anno_Stat{$_}{0}\t$Anno_Stat{$_}{1000}\t$Anno_Stat{$_}{2000}\t$Anno_Stat{$_}{3000}\t$Anno_Stat{$_}{6000}\n";
	}
	else {
		print STAT "All_Annotated\t$Anno_Stat{$_}{Anno}\t$Anno_Stat{$_}{0}\t$Anno_Stat{$_}{1000}\t$Anno_Stat{$_}{2000}\t$Anno_Stat{$_}{3000}\t$Anno_Stat{$_}{6000}\n";
	}
}

close STAT;


###### draw database annotation genes Venn Diagram ##############
&MKDIR ("$odir/Venn_list");
my @Venns_database = ("COG","KOG","Swissprot","NR","TrEMBL","KEGG");  
open LIST,">$odir/Venn_list/list" || die $!;
my $venn_index = 0;
for my $db (@Venns_database) {
	if (exists $Venn_list{$db} && $venn_index <=5) {
		$venn_index ++;
		open OUT,">$odir/Venn_list/$db" || die $!;
		print LIST "$odir/Venn_list/$db\n";
		foreach (keys %{$Venn_list{$db}}) {
			print OUT "$_\n";
		}
		close OUT;
	}
}
close LIST;

system "perl $Bin/Venn.pl -l $odir/Venn_list/list -od $odir -n $cds_key.Anno_Venn";
system "rm -r $odir/Venn_list";

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\n$programe_dir End Time :[$Time_End]\n\n";

###########subs
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

sub help
{
	print <<"	Usage End.";
	Description:
		Function : Use DEG_Analysis Out to Extract Annotation file and Draw Pictures(COG,Kegg..);
		Version  : $ver
		Writer   : mengf <mengf\@biomarker.com.cn>
		Usage    :
		-gene
		    Unigene fa file ,to state The gene length;
		-id
		    Annotation OUT DIR (eg. ./Anno_dir);
		-od
		    Anno Integerate and Stat dir (eg. /Anno_integrate);
		-key
		    mRNA fa name (mRNA seq file name)

	Usage End.

	exit;
}
