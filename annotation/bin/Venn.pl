#!/usr/bin/perl -w
use strict;

use lib '/home/mengfei/workdir/perl_module/lib/lib/perl5';
#use lib '/home/zhanghk/software/perl/lib/lib64/perl5';

use Statistics::R;
use Getopt::Long;
use File::Basename;
use Algorithm::Combinatorics qw(combinations permutations);

my ($list, $od, $pic, $help);

GetOptions(
	'l=s'  => \$list,
	'od=s' => \$od,
	'n=s' => \$pic,
	'h' => \$help
);

my $usage=<<INFO;
Usage:
    perl $0 [options]
Options:
	-l   <file>    contain group list filename, filename be used for name group
	-od  <dir>     output directory Venn diagram pdf
	-n   <file>    output file name, only support pdf format(.pdf add auto) 
	-h             show this help information
Example:
	perl $0 -l group.list -n rice
Author:
	zhanghaikuan\@berrygenomics.com
Note:
	Group number should not more than 5, otherwise use the first 5 for drawing
INFO

if($help || !$list || !$pic){
	die $usage;
}

open LIST,"<$list" or die "$!\n";
my @files = <LIST>;
chomp @files;
close LIST;

my %fills = (
	2=>"fill=c(\"blue\", \"red\")",
	3=>"fill=c(\"blue\", \"red\", \"green\")",
	4=>"fill=c(\"orange\", \"red\", \"green\", \"blue\")",
	5=>"fill=c(\"dodgerblue\", \"goldenrod1\", \"darkorange1\", \"seagreen3\", \"orchid3\")",
);

my %catcol = (
	2=>"cat.col=c(\"blue\", \"red\")",
	3=>"cat.col=c(\"blue\", \"red\", \"green\")",
	4=>"cat.col=c(\"orange\", \"red\", \"green\", \"blue\")",
	5=>"cat.col=c(\"dodgerblue\", \"goldenrod1\", \"darkorange1\", \"seagreen3\", \"orchid3\")",
);

my $gnumber = 0; #keep group number
my %hash = ();
my $category = "c(";

#############################start get draw area para##################################################
my $area = "";

for(my $i=0; $i<5 && $i<=$#files; $i++){
	open IN,"<$files[$i]" or die "$!\n";
	my $tile = basename $files[$i];
	$category.="\"$tile\",";
	while(<IN>){
		chomp;
		$hash{$i+1}->{$_} = 1;
	}
	$gnumber++;
}
chop($category);
$category.=")";

if($gnumber <=1){
	die "Sorry, only one group can't support!\n";
}

my @arry = 1..$gnumber;

foreach my $k (@arry){
	$area.=join("","area$k=",scalar(keys %{$hash{$k}}),",");
}

for(my $i=2; $i<=$gnumber; $i++){
	my $iter = combinations(\@arry, $i);
	while (my $c = $iter->next) {
		my $count = intersection($c, \%hash);
		$area.=join("", "n", @$c, "=$count,");
	}
}

$area=~s/n12=/cross.area=/ if($gnumber==2);

#############################end get draw area para##################################################

#############################start determin draw function##############################################

my %funs = (
	2=>"draw.pairwise.venn",
	3=>"draw.triple.venn",
	4=>"draw.quad.venn",
	5=>"draw.quintuple.venn",
);
my $draw_fun = $funs{$gnumber};

#############################end determin draw function################################################

chop($area);
my $R = Statistics::R->new();
$R->run(qq`library(VennDiagram)`);
$R->run(qq`png("$od/$pic.png")`);
$R->run(qq`venn.plot<-$draw_fun($area,$category,$fills{$gnumber},$catcol{$gnumber}, margin = 0.12,lty="blank", cex = 1.5, cat.cex=1.5)`);
$R->run(qq`tiff(filename = "$od/$pic.tiff", compression = "lzw")`);
$R->run(qq`grid.draw(venn.plot)`);
$R->run(qq`dev.off()`);
$R->stop();

sub intersection{
	my ($tags, $hs) = @_;
	my %temp = ();
	my $value = 0;
	foreach my $t (@$tags){
		foreach my $k (keys %{$hs->{$t}}){
			$temp{$k}+=1;
		}
	}
	my $len = scalar(@$tags);
	foreach my $k (keys %temp){
		if($temp{$k} == $len){
			$value++;
		}
	}
	return $value;
}
