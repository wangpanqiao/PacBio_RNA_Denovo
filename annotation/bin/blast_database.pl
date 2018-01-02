#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

my $program_name=basename($0);
my ($dir,$blastallP,$Evalue,$Nr,$Nt,$Swissprot,$TrEMBL,$Cog,$Kog,$Kegg,$Cpu,$Database_config,$Outdir);
my ($seq_file_name,$Verbose,$Help,$Queue);
GetOptions(
	"Q_dir:s"=>\$dir,
	"TAB:s"=>\$Outdir,
	"name:s"=>\$seq_file_name,
	"evalue:s"=>\$Evalue,
	"p:s"=>\$blastallP,
	"cfg:s"=>\$Database_config,
	"nr"=>\$Nr,
	"nt"=>\$Nt,
	"swissprot"=>\$Swissprot,
	"trembl"=>\$TrEMBL,
	"cog"=>\$Cog,
	"kog"=>\$Kog,
	"kegg"=>\$Kegg,
	"cpu:i"=>\$Cpu,
	"queue:s"=>\$Queue,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);

$Queue ||= "all.q";

###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n$program_name Start Time :[$Time_Start]\n\n";


#################

my %blast_config;
parse_config($Database_config,\%blast_config);

my $notename=`hostname`;chomp $notename;

my $blast_parser = "$Bin/blast_parser.pl";
my $nr_path = $blast_config{nr};
my $nt_path = $blast_config{nt};
my $swissprot_path = $blast_config{Swissprot};
my $trembl_path = $blast_config{TrEMBL};
my $cog_path = $blast_config{Cog};
my $kog_path = $blast_config{Kog};
my $kegg_path = $blast_config{Kegg};

my $blast_shell_file = "$dir/work_sh/$seq_file_name.blast.sh";
my @subfiles;
@subfiles = glob("$dir/$seq_file_name.div/*.fa");

############creat shell file
open OUT,">$blast_shell_file" || die "fail $blast_shell_file";
foreach my $subfile (@subfiles) {
	print OUT "/share/public/software/blast/blast-2.2.20/bin/blastall -b 100 -v 100 -p $blastallP -e $Evalue -F F -d $nr_path -i $subfile -m 7 -a 2 -o $subfile.NR.blast.xml && \n" if(defined $Nr);
	print OUT "/share/public/software/blast/blast-2.6.0/ncbi-blast-2.6.0+/bin/blastn  -evalue $Evalue -num_threads 16 -outfmt 0  -db $nt_path -query $subfile  -out $subfile.NT.blast && \n" if(defined $Nt);
	print OUT "/share/public/software/blast/blast-2.2.20/bin/blastall -b 100 -v 100 -p $blastallP -e $Evalue -F F -d $swissprot_path -i $subfile -m 7 -a 2 -o $subfile.Swissprot.blast.xml && \n" if(defined $Swissprot);
	print OUT "/share/public/software/blast/blast-2.2.20/bin/blastall -b 100 -v 100 -p $blastallP -e $Evalue -F F -d $trembl_path -i $subfile -a 2 -o $subfile.TrEMBL.blast && \n" if(defined $TrEMBL);
	print OUT "/share/public/software/blast/blast-2.2.20/bin/blastall -b 100 -v 100 -p $blastallP -e $Evalue -F F -d $kegg_path -i $subfile -m 9 -a 2 -o $subfile.KEGG.blast.tab && \n" if(defined $Kegg);
	print OUT "/share/public/software/blast/blast-2.2.20/bin/blastall -b 100 -v 100 -p $blastallP -e $Evalue -F F -d $cog_path -i $subfile -a 2 -o $subfile.COG.blast && \n" if(defined $Cog);
	print OUT "/share/public/software/blast/blast-2.2.20/bin/blastall -b 100 -v 100 -p $blastallP -e $Evalue -F F -d $kog_path -i $subfile -a 2 -o $subfile.KOG.blast && \n" if(defined $Kog);
}
close OUT;

####################run the shell file
&Cut_shell_qsub("$blast_shell_file",$Cpu,"10G","$Queue");

####################

################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\n$program_name End Time :[$Time_End]\n\n";

####################################################
################### Sub Routines ###################
####################################################

##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $blast_config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if(/^$/ or /^\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$blast_config_p->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}

	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime
{ # &Runtime($BEGIN);
        my ($t1)=@_;
        my $t=time()-$t1;
        print "\nBlast Total elapsed time: ${t}s\n";
}

sub Show_time{#Show time
	my($step)=@_;
	my $time = sub_format_datetime(localtime(time()));
	print "\n$step Time :[$time]\n\n";
}

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
		@div;
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
