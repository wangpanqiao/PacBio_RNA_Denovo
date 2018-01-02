#/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Term::ANSIColor qw(:constants);

my($in, $outdir, $run, $help);

GetOptions(
	"i|in=s" => \$in,
	"o|outdir=s" => \$outdir,
	"r|run" => \$run,
	"h|help" => \$help
);

$outdir ||= ".";
my $formatdb = "/share/public/software/blast/blast-2.2.20/bin";
my $blast = "/share/public/software/ncbi-blast-2.2.25+/bin";
system("mkdir $outdir") unless (-d $outdir);
system("mkdir $outdir/result") unless (-d "$outdir/result");
&usage if(!$in || $help);

open SH, ">$outdir/AS.sh" or die "cannot open:$!";
print SH "#step 1: all vs all blast transcripts\n";
print SH "$formatdb/formatdb -i $in -p F\n";
print SH "$blast/blastn -query $in -db $in -out $outdir/all_vs_all_blast.xml -evalue 1e-10 -num_threads 14 -perc_identity 99 -outfmt 5\n";
print SH "#step 2: identity the AS event\n";
print SH "python $Bin/bin/1step.py $outdir/all_vs_all_blast.xml $outdir/temp1.xls\n";
print SH "python $Bin/bin/2.1step.py $outdir/temp1.xls $outdir/temp2.1.xls\n";
print SH "python $Bin/bin/2.2step.py $outdir/temp2.1.xls $outdir/temp2.2.xls\n";
print SH "python $Bin/bin/3step.py $outdir/temp2.2.xls $outdir/temp3.xls\n";
print SH "python $Bin/bin/4step.py $outdir/temp3.xls $outdir/temp4.xls\n";
print SH "python $Bin/bin/5step.py $outdir/temp4.xls $outdir/temp5.xls\n";
print SH "python $Bin/bin/6.1step.py $outdir/temp5.xls $outdir/temp6.1.xls\n";
print SH "python $Bin/bin/6.2step.py $outdir/temp6.1.xls $outdir/temp6.2.xls\n";
print SH "sed 's/No definition line//g' $outdir/temp6.2.xls > $outdir/temp6.3.xls\n";
print SH "sort -k 2,2V -k 3,3rn $outdir/temp6.3.xls > $outdir/temp6.4.xls\n";
print SH "awk '{if(! a[\$2]){print; a[\$2]++}}' $outdir/temp6.4.xls > $outdir/temp6.5.xls\n";
print SH "sort -k 1,1V -k 3,3rn $outdir/temp6.5.xls > $outdir/temp6.6.xls \n";
print SH "awk '{if(! a[\$1]){print; a[\$1]++}}' $outdir/temp6.6.xls > $outdir/temp7.xls\n";
print SH "#step 3: plot and summary the as results\n";
print SH "Rscript $Bin/bin/plot_as.R -i $outdir/temp7.xls -o $outdir/result\n";
close SH;

system("cd $outdir && qsub -cwd -l vf=10g $outdir/AS.sh") if($run);

sub usage{
	my $info =<<INFO;

	Usage:
		--i|in :           The transcript sequences for AS analysis.
		--o|outdir:        The outdir of results.
		--h|help:	   If set it, it will show the help message.
INFO
	print RED $info;
	exit;
}
