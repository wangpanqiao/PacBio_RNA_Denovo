#!/usr/bin/perl -w
use strict;
use FindBin '$Bin';
use Getopt::Long;

my $usage =
"
Usage:

	Options:
	-fasta     the fasta files
        -outdir    the outdir of files
        -name      file name
	-h     Help

For example:
	 perl -f gene.fa -o ./
";

my ($fasta, $outdir, $name, $help);
GetOptions(
	"h|help"=>\$help,
        "o|outdir=s"=>\$outdir,
        "n|name=s"=>\$name,
	"f|fasta=s"=>\$fasta,
);
if(!$fasta || !$outdir || !$name || $help){
	die "$usage\n";
}
my $primer3 = '/share/public/software/primer3/primer3_0_6/src/primer3_core';
system("cp $Bin/misa.ini $outdir/") unless(-e "$outdir/misa.ini");
system("cp $fasta $outdir/$name");
system("perl $Bin/name_line_convert.pl $outdir/$name");
system("mv $outdir/$name.single.fa $outdir/$name");
system("perl $Bin/misa.pl $outdir/$name");
system("perl $Bin/p3_in.pl $outdir/$name.misa");
system("$primer3 -default_version=1 -output=$outdir/$name.p3out $outdir/$name.p3in");
system("perl $Bin/p3_out.pl $outdir/$name.p3out $outdir/$name.misa");
system("mv $outdir/$name.results $outdir/$name.results.xls");
system("perl $Bin/ssr_stat.pl $outdir/$name $outdir/$name.results.xls ssr_density.txt");
