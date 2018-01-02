#!/usr/bin/perl -w                                                                                                                             
use strict;
use Term::ANSIColor qw(:constants);
	$Term::ANSIColor::AUTORESET=1;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

use FindBin qw($Bin);
use lib "$Bin/../lib/";
use FASTA;
use TOOLS;
use MAKEFILE;

my ($sample_info, $outdir, $xml, $species, $help);

GetOptions( "sample|s:s"  => \$sample_info,
			"xml|x:s"     => \$xml,
			"od|o:s"      => \$outdir,
			"species:s"   => \$species,
			"help|h|?"    => \$help
);

&usage if(!$sample_info || !$outdir || !$xml || $help);
$xml = abs_path($xml);
mkdir($outdir) unless(-d $outdir);
$species ||= "species";
my $config_ini  = "$outdir/assemble.ini";
-e $config_ini or `cp  $Bin/assemble.ini  $config_ini`;
my $config  = TOOLS::read_ini ($config_ini);
$$config{BIN} = $Bin;
my $samples_fofn = create_fofn($sample_info, $outdir);

my @make;
my $j = 0;

foreach my $sample(keys %$samples_fofn){
	$make[$j]{I} = "";
	$make[$j]{D} = "$sample";
	$make[$j]{O} = "$species\_$sample.fasta";
	#print "$samples_fofn->{$sample}\n";
	#print "$outdir/$sample\n";
	#print "$sample\n";
	#$make[$j]{C} = TOOLS::cmd ($config, 'SMRT_CMD', $samples_fofn->{$sample}, $outdir/$sample, $xml, $species, $sample);
	$make[$j]{C} = TOOLS::cmd ($config, 'SMRT_CMD', $samples_fofn->{$sample}, "$outdir/$sample", $xml, $species, $sample);

	$j++;
}

MAKEFILE::write_makefile ($outdir, 'assemble', \@make);

sub create_fofn{
	my $file = shift;
	my $dir  = shift;
	my %hash;
	open IN, "<$file" or die "[ERROR] Sample information file open failed: $!\n";
	while(<IN>){
		chomp;
		next if(/^#/ || /^\s*$/);
		my @eles = split /[\s,]+/, $_;
		for(my $i=1; $i<=$#eles; $i++){
			push @{$hash{$eles[0]}}, $eles[$i];
		}
	}
	close IN;
	my %sample_fofn;
	foreach my $sample (keys %hash){
		mkdir("$dir/$sample") unless(-d "$dir/$sample");
		my %h5;
		foreach my $d (@{$hash{$sample}}){
			my @h5s = `find $d -name "*bax.h5"`;
			chomp @h5s;
			map{$h5{$_}=1} @h5s;
		}
		open FOFN, ">$dir/$sample/$sample.fofn" or die "[ERROR] create fofn file failed: $!\n";
		print FOFN join("\n", keys %h5)."\n";
		close FOFN;
		$sample_fofn{$sample} = "$dir/$sample/$sample.fofn";
	}
	return \%sample_fofn;
}

sub usage{
	my $info=<<INFO;

Usage   : perl $0 [options]
Options : s|sample <s>    Sample information(sample_name  h5_dir1,h5_dir2)
          species  <s>    Assembly species name, default(species)
          x|xml    <s>    Parameters for SMRT in XML format
          o|od     <s>    Assembly result output directory
          h|help|?        Show this help information
Example: perl $0 -sample -xml -od 

INFO
	die GREEN $info;
}
