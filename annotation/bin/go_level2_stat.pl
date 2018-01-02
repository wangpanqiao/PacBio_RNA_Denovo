#!/usr/bin/perl -w
#use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0";
my @Times = localtime();
#######################################################################################
my $Time_Start = sub_format_datetime(localtime(time()));
print STDOUT "$Script start at:[$Time_Start]\n";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($anno,$lable,$foutdir,$fkey,$nodraw,$note,$partdraw,$color_index);
GetOptions(
				"help|h|?" =>\&USAGE,
				"gene2GO:s"=>\$anno,
				"lable:s"=>\$lable,
				"od:s"=>\$foutdir,
				"prefix:s"=>\$fkey,
				"note:s"=>\$note,
				"partdraw"=>\$partdraw,
				"color:i"=>\$color_index,
				) or &USAGE;
&USAGE unless ($anno and $foutdir and $fkey and $lable);

######## paras input ##############
my @annos = split /,/,$anno;
for (0..$#annos) {
	$annos[$_] = Cwd::realpath($annos[$_]);
}

my @lables = split /,/,$lable;

if ($#annos != $#lables) {
	print "input gene2GO files unequal to out prefix or graph lables \n\n";
	die;
}
if ($#annos > 3) {
	print "the program only support 4 gene2GO files graph\n\n";
	die;
}

$foutdir = Cwd::realpath($foutdir);
mkdir $foutdir if (!-d $foutdir);
mkdir "$foutdir/tmp" if (!-d "$foutdir/tmp");
my $nowdir = Cwd::realpath();
$color_index = $color_index || 1;

######## database file############
my %three_main_ontology = (
		"GO:0008150" => "biological_process",
		"GO:0005575" => "cellular_component",
		"GO:0003674" => "molecular_function",
		);

my %part_mark = (
		"biological_process" => "BP",
		"cellular_component" => "CC",
		"molecular_function" => "MF",
		);

my %level2term_files = (
		"biological_process" => "/share/work2/liuying/condatabase/GO/BP_level2.terms",
		"cellular_component" => "/share/work2/liuying/condatabase/GO/CC_level2.terms",
		"molecular_function" => "/share/work2/liuying/condatabase/GO/MF_level2.terms",
		);

my %level2terms;
&load_level2_ (\%level2term_files,\%level2terms);

######## stat and graph list #########
my %all_level2stat;
my %terms_occur;
my %anno_sum;
my %three_main_ontology_stats;

for my $index (0..$#annos) {

	my $sum = `wc -l $annos[$index]`;chomp $sum;
	$sum =~ s/^\s+//;
	$sum =~ s/\s.*$//;
	$anno_sum{$lables[$index]} = $sum;
	######## stat GO terms ####
	system "cd $foutdir/tmp && Rscript $Bin/go_terms_stat.r $annos[$index] $lables[$index]";

	######## extract level2 stats ####
	my %term_stats = (
			"biological_process" => "$foutdir/tmp/$lables[$index]"."_BP_TermStat.xls",
			"cellular_component" => "$foutdir/tmp/$lables[$index]"."_CC_TermStat.xls",
			"molecular_function" => "$foutdir/tmp/$lables[$index]"."_MF_TermStat.xls",
			);

	my %level2terms_stats;
	&extract_stat_ (\%term_stats,\%three_main_ontology,\%level2terms,\%three_main_ontology_stats,\%level2terms_stats);

	######## stat results ####
	foreach my $onto (sort keys %level2terms_stats) {
		foreach my $term (keys %{$level2terms_stats{$onto}}) {
			$terms_occur{$onto}{$term} = $level2terms{$onto}{$term};
			$all_level2stat{$lables[$index]}{$onto}{$term} = $level2terms_stats{$onto}{$term};
		}
	}
}

######## output stat ###############
my %percent;
my %ontology_stat;

open OUT,">$foutdir/$fkey"."_GOTermStat_level2.xls" || die $!;
print OUT "GOTerm\tGO_def\t";
print OUT join "\t",@lables,"\n";
foreach my $onto (sort keys %terms_occur) {
	print OUT "$onto\t$three_main_ontology_stats{$onto}\n";
	foreach my $term (sort keys %{$terms_occur{$onto}}) {
		$ontology_stat{$onto} ++;

		my $str = "$term\t$terms_occur{$onto}{$term}\t";
		for my $index (0..$#lables) {
			if (exists $all_level2stat{$lables[$index]}{$onto}{$term}) {
				$str .= "$all_level2stat{$lables[$index]}{$onto}{$term}\t";
				$percent{$onto}{$term}{$lables[$index]} = sprintf "%.4f",($all_level2stat{$lables[$index]}{$onto}{$term}/$anno_sum{$lables[$index]}+0.001)*100;
			}
			else {
				$str .= "0\t";
				$percent{$onto}{$term}{$lables[$index]} = "0.1";
			}
		}
		$str =~ s/\s+$//;
		print OUT "$str\n";
	}
}
close OUT;

######## output graph list file ############
my $dissvg_pl="$Bin/distributing_svg_4.74/distributing_svg.pl";
my $svgxxx="$Bin/distributing_svg_4.74/svg2xxx_release/svg2xxx";

my $note_type = "GO Classification";
$note_type = "$note GO Classification" if (defined $note);
my $graph_move=0.1;

if (!defined $nodraw) {
	my %colors = (
			"1" => ["#0067A6","#F2572D","#00ABD8","#EFC028"],
			"2" => ["#FF6666","#3399CC","#FF9900","#009966"],
			"3" => ["#990033","#CCFF66","#FF9900","#333399"],
			);

	my @rect_width=("0.55","0.3","0.25","0.2");
	my @moveper=("0.4","0.3","0.125","0.1");

	#--------------------------------------------
	######## prapare for graph infos ########
	my $x_end = 0;
	my $x_scale_axis;
	my %lables_info;
	my %part_lables_info;
	my %part_x_end;
	my %part_x_scale_axis;

	foreach my $onto (sort keys %percent) {
		my $part_x = 0;
		my $part_x_scale;
		foreach my $term (sort keys %{$percent{$onto}}) {
			$x_scale_axis .= "$terms_occur{$onto}{$term}\n";
			$part_x_scale .= "$terms_occur{$onto}{$term}\n";
			for my $index (0..$#lables) {
				$lables_info{$lables[$index]} .= "$x_end:$percent{$onto}{$term}{$lables[$index]}\n";
				$part_lables_info{$lables[$index]}{$onto} .= "$part_x:$percent{$onto}{$term}{$lables[$index]}\n";
			}
			$x_end ++;
			$part_x ++;
		}
		$part_x_end{$onto} = $part_x;
		$part_x_scale_axis{$onto} = $part_x_scale;
	}

	my $x_scale;
	foreach my $onto (sort keys %ontology_stat) {
		$x_scale .= "$ontology_stat{$onto}:$onto\n";
	}

	#--------------------------------------------

	my $graph_svg_info;
	$graph_svg_info = &graph_out($x_end, $rect_width[scalar(@lables)-1], $moveper[scalar(@lables)-1]);
	$graph_svg_info .= $x_scale_axis;
	$graph_svg_info .= ":END\n";

	$graph_svg_info .= "Group:\n"."$x_scale".":End\n\n";

	for my $index (0..$#lables) {
		$graph_svg_info .= "Color:${$colors{$color_index}}[$index]\n";
		$graph_svg_info .= "YMark:r\n";
		$graph_svg_info .= "Start:0\n";
		$graph_svg_info .= "End:3\n";
		$graph_svg_info .= "Step:1\n";

		$graph_svg_info .= "Scale:\n";
		$graph_svg_info .= &sub_def_Scale_ ($anno_sum{$lables[$index]});
		$graph_svg_info .= ":End\n";

		$graph_svg_info .= "Mark:$lables[$index]\n";
		$graph_svg_info .= "$lables_info{$lables[$index]}\n\n";
	}

	my $graph_prefix = "$fkey"."_GOTermStat_level2";
	open OUT,">$foutdir/$graph_prefix.svg.list" || die $!;
	print OUT "$graph_svg_info";
	close OUT;

	system "perl $dissvg_pl $foutdir/$graph_prefix.svg.list $foutdir/$graph_prefix.svg";
	&tuning_svg("$foutdir/$graph_prefix.svg");
	system "mv $foutdir/$graph_prefix.svg.raw $foutdir/$graph_prefix.svg";
	chdir($foutdir);
	system "perl $svgxxx $graph_prefix.svg -t png";
	system "perl $svgxxx $graph_prefix.svg -t pdf";
	chdir($nowdir);


	if (defined $partdraw) {
		foreach my $onto (sort keys %ontology_stat) {
			my $graph_part_svg;
			$graph_part_svg = &graph_part_out ($part_x_end{$onto}, $rect_width[scalar(@lables)-1], $moveper[scalar(@lables)-1], $note_type);
			$graph_part_svg .= "$part_x_scale_axis{$onto}";
			$graph_part_svg .= ":END\n";

			$graph_part_svg .= "Group:\n"."$ontology_stat{$onto}:$onto\n".":END\n\n";

			for my $index (0..$#lables) {
				$graph_part_svg .= "Color:${$colors{$color_index}}[$index]\n";
				$graph_part_svg .= "YMark:r\n";
				$graph_part_svg .= "Start:0\n";
				$graph_part_svg .= "End:3\n";
				$graph_part_svg .= "Step:1\n";

				$graph_part_svg .= "Scale:\n";
				$graph_part_svg .= &sub_def_Scale_ ($anno_sum{$lables[$index]});
				$graph_part_svg .= ":End\n";

				$graph_part_svg .= "Mark:$lables[$index]\n";
				$graph_part_svg .= "$part_lables_info{$lables[$index]}{$onto}\n\n";
			}

			my $praph_part_prefix = "$fkey"."_$part_mark{$onto}"."_GOTermStat_level2";
			open OUT,">$foutdir/$praph_part_prefix.svg.list" || die $!;
			print OUT "$graph_part_svg";
			close OUT;

			system "perl $dissvg_pl $foutdir/$praph_part_prefix.svg.list $foutdir/$praph_part_prefix.svg";
			&tuning_svg("$foutdir/$praph_part_prefix.svg");
			system "mv $foutdir/$praph_part_prefix.svg.raw $foutdir/$praph_part_prefix.svg";
			chdir($foutdir);
			system "perl $svgxxx $praph_part_prefix.svg -t pdf";
			chdir($nowdir);
		}
	}
}

#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub load_level2_ {
	my ($file,$class) = @_;

	my %files = %$file;
	foreach my $ontology (sort keys %files) {
		open IN,"$files{$ontology}" || die $!;
		while (<IN>) {
			chomp;
			s/^\s+//;s/\s+$//;
			next if ($_ !~ /^GO:/);
			my ($term,$def) = split /\t+/,$_;
			$class -> {$ontology} {$term} = $def;
		}
		close IN;
	}
}

sub extract_stat_ {
	my ($file,$level1term,$level2term,$level1termstat,$level2termstat) = @_;

	my %files = %$file;
	foreach my $ontology (sort keys %files) {
		open IN,"$files{$ontology}" || die $!;
		while (<IN>) {
			chomp;
			s/^\s+//;s/\s+$//;
			next if ($_ !~ /^GO:/);
			my ($goterm,$anno) = split /\t+/,$_;
			if (exists $level1term -> {$goterm} ) {
				$level1termstat -> {$ontology} = $anno;
			}
			if (exists $level2term -> {$ontology}{$goterm} ) {
				$level2termstat -> {$ontology}{$goterm} = $anno;
			}
		}
		close IN;
	}
}

sub sub_def_Scale_ {
	my ($limit) = @_;

	my $out;
	my $len = length($limit);
	my $f=0;
	my $wei;
	for (my $j=3;$j>=0;$j--){
		$wei = $len-$j;
		if ($wei<0)
		{
			$out .= "0\n";
			$f=1;
		}
		elsif($wei==0)
		{
			if ($f==1) 
			{
				$out .= "1\n";
			}
			else
			{
				$out .= "0\n";
			}
		}
		else 
		{
			$out .= substr($limit,0,$wei);
			$out .= "\n";
		}
	}
	return ($out);
}

sub graph_out {
	my ($x_end,$width,$moveper)=@_;
	$moveper+=$graph_move;
	my $dis =<<"OUT";
Type:Simple
Width:3000
Height:800
BothYAxis:1
MultiRY:1
ScaleLen:8
WholeScale:1
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
YStep:33.3
RYStep:6
XScalePos:0.5
XScaleRoate:80
XStart:0
YStart:0.1
RYStart:0.1
XEnd:$x_end
YEnd:100
RYEnd:100
YNeedLog:10
MarkNoBorder:1
WholeScale:0.95
Note:$note_type 
X:
Y:Percent of genes
RY:Number of genes
XUnit:1
Scale:
OUT

return $dis;
}

sub graph_part_out ($$$$) {
	my ($x_end,$width,$moveper,$part_note)=@_;
	$moveper+=$graph_move;
	my $part_width = 150 * ($x_end - 1);
	my $dis =<<"OUT";	
Type:Simple
Width:$part_width
Height:800
BothYAxis:1
MultiRY:1
ScaleLen:8
WholeScale:1
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
YStep:33.3
RYStep:6
XScalePos:0.5
XScaleRoate:80
XStart:0
YStart:0.1
RYStart:0.1
XEnd:$x_end
YEnd:100
RYEnd:100
YNeedLog:10
MarkNoBorder:1
WholeScale:0.95
Note: $part_note
X:
Y:Percent of genes
RY:Number of genes
XUnit:1
Scale:
OUT

return $dis;
}

sub tuning_svg {
	my $raw=shift;
	open (IN,$raw) ||die;
	open (OUT,">$raw.raw")||die;
	while(<IN>)
	{
		chomp;
		if(/^\<svg width=\"([.\d]+)\"/){
			s/($1)/$1+150/e;
		}
		print OUT "$_\n";
	}
	close IN;
	close OUT;
}

sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {
	my $usage=<<"USAGE";
	Program: $Script
	Version: $version: [2015/1/22]
	Contact: Mengfei <fred_routine\@163.com.cn>
	Description: Extract KEGG Pathway and KOs file From kobas anno file;

	Usage:
		-gene2GO     The blast_tab Out of seq with KEGG                     must be given
		             support maximum 4 files graph, seperated by comma;

		-lable       graph lable of go histogram                            must be given
		             is to distinguish the data in the plots in different color.
		             [the number of lables must equal that of input files, seperated by comma;]

		-od          Output dir                                             must be given

		-prefix      Prefix of OUT files (eg. prefix.BP_level2.stat.xls)    must be given

		-note        graph note on top of go graph                          optical

		-partdraw    graph three ontologys partly                           optical

		-color       hist color sets scheme (1 or 2 or 3)                   optical
		             default 1 ("#CC0033","#006699","#FF9933","#663399");

		-h          Help document
USAGE
	print $usage;
	exit;
}
