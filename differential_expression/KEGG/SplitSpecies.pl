#! /urs/bin/perl -w 
#use strict;
use FindBin qw($Bin $Script);
my %hash;
my %con;
open ALL,"$ARGV[0]";
while (<ALL>){
	chomp;
	next if /^#/;
	my @tem;
	my $spe;
	@tem=split /\t/;
	$spe=(split /:/,$tem[1])[0];
	if (exists $hash{$spe}){
		$con{$spe}="$con{$spe}$_\n";
	}else{ 
		$con{$spe}="$_\n";
		$hash{$spe}=1;
		}
}
close ALL;

open SH, "> $ARGV[1]/anno.sh";
foreach my $spe (sort keys %hash){
	open $spe,"> $ARGV[1]/$spe\_blast.txt";
	print $spe "$con{$spe}";
	print SH "export PYTHONPATH=/share/work3/minghb/software/kobas-3.0/src:\$PYTHONPATH && python $Bin/annotate.py -i $ARGV[1]/$spe\_blast.txt -t blastout:tab -s $spe -o $ARGV[1]/$spe.annotate && rm $ARGV[1]/$spe\_blast.txt &&\n";   
	close $spe;
	}
close SH;

system "ssh login-0-0  perl $Bin/qsub-sge.pl --queue all.q  --resource vf=5.0G $ARGV[1]/anno.sh";
