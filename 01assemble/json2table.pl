#!/usr/bin/perl -w
use strict;
use JSON;
use Encode;
use File::Basename;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my ($assemble_dir, $outdir, $help);

GetOptions( "ad:s" => \$assemble_dir,
			"od:s" => \$outdir,
			"h|help|?" => \$help,
);

&usage if(!$assemble_dir || !$outdir || $help);

my @report_filter = `find $assemble_dir -name filter_reports_filter_stats.json`;
my @report_polish = `find $assemble_dir -name polished_report.json`;
chomp @report_filter;
chomp @report_polish;

my $head1;
my @filter_lines;
foreach my $file (@report_filter){
	$file =~ /assemble\/([^\/]+)/;
	my $name = $1;
	my ($temp, $line) = filter($name, $file);
	$head1 = $temp;
	push @filter_lines, $line;
}
my $head2;
my @polish_lines;
foreach my $file (@report_polish){
	$file =~ /assemble\/([^\/]+)/;
	my $name = $1;
	my ($temp, $line) = polish($name, $file);
	$head2 = $temp;
	push @polish_lines, $line;
}

outfile($head1, \@filter_lines, "$outdir/filter_stat.txt");
outfile($head2, \@polish_lines, "$outdir/polish_stat.txt");

## Sub Functions
#
sub filter{
	my ($sample_name, $json_file) = @_;
	local $/;

	open FH, "<$json_file";
	my $json_text = <FH>;
	close FH;

	my $perl_scalar = from_json( $json_text, { utf8  => 1 } );

	my $p1 = ${${$perl_scalar->{'tables'}}[0]->{'columns'}}[0]->{'values'};
	my $v1 = ${${$perl_scalar->{'tables'}}[0]->{'columns'}}[1]->{'values'};
	my $v2 = ${${$perl_scalar->{'tables'}}[0]->{'columns'}}[2]->{'values'};
	my $p2 = $perl_scalar->{'attributes'};

	my @head;
	foreach my $value (@$p1){
		push @head, $value;
	}

	my @content;
	my $flag = 0;
	foreach my $value (@$v1){
		push @content, "$value/$$v2[$flag]";
		$flag++;
	}

	foreach my $atr (@$p2){
		push @head, $atr->{'name'};
		push @content, "-/$atr->{'value'}";
	}
	splice(@head, 5);
	splice(@content, 5);
	return (join("\t", "SampleName", @head), join("\t", $sample_name, @content));
}

sub polish{
	my ($sample_name, $json_file) = @_;
	local $/;

	open FH, "<$json_file";
	my $json_text = <FH>;
	close FH;

	my $perl_scalar = from_json( $json_text, { utf8  => 1 } );

	my @head;
	my @content;
	my $p1 = $perl_scalar->{'attributes'};
	foreach my $atr (@$p1){
		push @head, $atr->{'name'};
		push @content, "$atr->{'value'}";
	}
	return (join("\t", "SampleName", @head), join("\t", $sample_name, @content));

}

sub outfile{
	my ($head, $lines, $fname) = @_;
	open OUT, ">$fname" or die "$!\n";
	#print "$head\n";die;
	print OUT $head."\n";
	foreach my $l (@$lines){
		print OUT $l."\n";
	}
	close OUT;
}

sub usage{
	my $info =<<INFO;
Usage   : perl $0 [options]
Options : ad  <dir>   Assembly result directory of all samples
          od  <dir>   Output statistic result directory
          h|help|?    Show this help information
Example : perl $0 -ad -od
INFO
	die $info;
}

#
#my $json = JSON->new->allow_nonref;
#my $pretty_printed;
#$pretty_printed = $json->pretty->encode( $perl_scalar );
#
#print $pretty_printed;


