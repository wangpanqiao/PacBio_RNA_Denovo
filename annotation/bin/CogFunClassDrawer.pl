#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="2.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$db,$PNG,$note);
GetOptions(
				"help|?"  => \&USAGE,
				"i:s" => \$fIn,
				"note:s" => \$note,
				"db:s" => \$db,
				"o:s" => \$fOut,
				"png"=>\$PNG,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $db);

my $ConvertSVGtoPNG="$Bin/distributing_svg_4.74/svg2xxx_release/svg2xxx";

# ------------------------------------------------------------------
# Class Anno
# ------------------------------------------------------------------
my %Category = (
	"Information Storage and Processing" => ["J","A","K","L","B"],
	"Cellular Processes and Signaling" => ["D","Y","V","T","M","N","Z","W","U","O"],
	"Metabolism" => ["C","G","E","F","H","I","P","Q"],
	"Poorly Characterized" => ["R","S"],
);

my %Class = (
	"J" => "Translation, ribosomal structure and biogenesis",
	"A" => "RNA processing and modification",
	"K" => "Transcription",
	"L" => "Replication, recombination and repair",
	"B" => "Chromatin structure and dynamics",
	"D" => "Cell cycle control, cell division, chromosome partitioning",
	"Y" => "Nuclear structure",
	"V" => "Defense mechanisms",
	"T" => "Signal transduction mechanisms",
	"M" => "Cell wall/membrane/envelope biogenesis",
	"N" => "Cell motility",
	"Z" => "Cytoskeleton",
	"W" => "Extracellular structures",
	"U" => "Intracellular trafficking, secretion, and vesicular transport",
	"O" => "Posttranslational modification, protein turnover, chaperones",
	"C" => "Energy production and conversion",
	"G" => "Carbohydrate transport and metabolism",
	"E" => "Amino acid transport and metabolism",
	"F" => "Nucleotide transport and metabolism",
	"H" => "Coenzyme transport and metabolism",
	"I" => "Lipid transport and metabolism",
	"P" => "Inorganic ion transport and metabolism",
	"Q" => "Secondary metabolites biosynthesis, transport and catabolism",
	"R" => "General function prediction only",
	"S" => "Function unknown",
);

my $TILE="$db Function Classification of Consensus Sequence";
$TILE="$db Function Classification of $note Unigenes" if (defined $note);
my $Xtile="Function Class";
my $Ytile="Frequency";

my %colors = (
	"Information Storage and Processing" => "#CD3333",
	"Cellular Processes and Signaling" => "#4876FF",
	"Metabolism" => "#FFD700",
	"Poorly Characterized" => "#00CD66",
);

# ------------------------------------------------------------------
# Get Data
# ------------------------------------------------------------------
my %Data;
&LoadData($fIn,\%Data);
open (OUT,">$fIn.stat.xls") || die "can't open file  $fIn.stat\n ";
&OutHash(\%Data);
close(OUT);

# ------------------------------------------------------------------
# SVG parameters
# ------------------------------------------------------------------
my $LeftMargin=20;
my $RightMargin=20;
my $TopMargin=20;
my $DownMargin=40;

my $MarkMargin = 60;
my $MarkFontSize = 18;

my $TileFontSize=25;
my $TileHeight=20;
my $TileLeftSpace=150;
my $TileFontColor="black";
my $TileDownSpace=20;
my $TileFont="ArialNarrow-Bold";

#my $LegendFont="Times New Roman";
my $LegendFont="ArialNarrow-Bold";
my $LegendFontSize=15;
my $LegendFontColar="black";
my $LegendIntervalSpace=20;
my $LegendLeftspace=20;
my $LegendWidth=850;

my $BarWidth=10;
my $BarInterval=5;

my $BoxWidth=($BarInterval+$BarWidth)*(scalar keys %Class)+$BarInterval;
my $BoxHeight=0;  #will be calculated later.

my $ScaleNum=6;
my $XYScaleLineLength=5;
my $XYScaleFontWidth=30;
my $XYScaleFontSize=15;
my $XYTileHeight=45;
my $XYTileFontSize=20;


my $PaperHeight=$TopMargin+$DownMargin+$TileHeight+$TileDownSpace+$MarkMargin+$LegendIntervalSpace*(scalar keys %Class);
my $PaperWidth=$LeftMargin+$RightMargin+$LegendWidth+$LegendLeftspace+$XYScaleLineLength+$XYScaleFontWidth+$XYTileHeight;

$BoxHeight=$PaperHeight-$TopMargin-$DownMargin-$TileHeight-$XYTileHeight-$TileDownSpace-$MarkMargin;

print "$PaperWidth,$BoxHeight\n";

# ------------------------------------------------------------------
# Set Scale base on data
# ------------------------------------------------------------------
my ($Max)=sort {$b <=> $a} values %Data;
my $ScaleReal=int(int($Max/0.9)/($ScaleNum));
my $Scale=(int($ScaleReal/10**(length($ScaleReal)-1))+1)*10**(length($ScaleReal)-1);

# ------------------------------------------------------------------
# Drawing
# ------------------------------------------------------------------
open (SVG,">",$fOut) or die $!;
print SVG &svg_paper($PaperWidth,$PaperHeight),"\n";

print SVG &svg_txt_litc($LeftMargin+$TileLeftSpace,$TopMargin+$TileHeight,$TileFontSize,$TileFontColor,$TileFont,$TILE);
print SVG &svg_rect_litc($LeftMargin+$XYScaleLineLength+$XYScaleFontWidth+$XYTileHeight,$TopMargin+$TileHeight+$TileDownSpace,$BoxWidth,$BoxHeight,"white","black",1);

my $curr=0;
my $offset=10;
foreach my $cate (sort keys %Category) {

	foreach my $ClassAnno (@{$Category{$cate}}) {

		my $value=exists $Data{$ClassAnno}?$Data{$ClassAnno}:0;
		$value=$value*$BoxHeight/($Scale*$ScaleNum);

		print SVG &svg_txt_litc($LeftMargin+$XYTileHeight+$XYScaleFontWidth+$XYScaleLineLength+$BoxWidth+$LegendLeftspace,$TopMargin+$TileHeight+$TileDownSpace+$LegendIntervalSpace*($curr+1)-$LegendFontSize/2,$LegendFontSize,$LegendFontColar,$LegendFont,$ClassAnno.":  ".$Class{$ClassAnno});
		
		print SVG &svg_txt_litc($LeftMargin+$XYTileHeight+$XYScaleFontWidth+$XYScaleLineLength+$BarInterval+($BarInterval+$BarWidth)*$curr,$offset+$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight+5,$LegendFontSize,"black",$LegendFont,$ClassAnno);

		print SVG &svg_rect_litc($LeftMargin+$XYTileHeight+$XYScaleFontWidth+$XYScaleLineLength+$BarInterval+($BarInterval+$BarWidth)*$curr,$TopMargin+$TileHeight+$TileDownSpace+($BoxHeight-$value),$BarWidth,$value,$colors{$cate},"black",1);

		$curr++;
	}
}

for (my $i=0;$i<=$ScaleNum ;$i++) {
	print SVG &svg_line($LeftMargin+$XYTileHeight+$XYScaleFontWidth+$XYScaleLineLength,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight-$i*$Scale*($BoxHeight/($Scale*$ScaleNum)),$LeftMargin+$XYTileHeight+$XYScaleFontWidth,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight-$i*$Scale*($BoxHeight/($Scale*$ScaleNum)),"black",1);

	print SVG &svg_txt_litc($LeftMargin+$XYTileHeight,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight-$i*$Scale*($BoxHeight/($Scale*$ScaleNum))+5,$XYScaleFontSize,"black",$LegendFont,$Scale*$i);

}

my @marks = sort keys %Category;
my $MarkXStep = $PaperWidth/2;
my $MarkYStep = ($MarkMargin)/(scalar(@marks));
my $SpaceMarkandXTile = $DownMargin/2;
print SVG &svg_rect_litc($PaperWidth*1/8-10,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight+$XYTileHeight+$SpaceMarkandXTile+$MarkYStep-10,$PaperWidth*8/12,$MarkMargin,"white","black",2);
for (my $i=0;$i<=$#marks ;$i+=2) {
	print SVG &svg_rect_litc($PaperWidth*1/8,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight+$XYTileHeight+$SpaceMarkandXTile+$MarkYStep*($i+1),10,10,$colors{$marks[$i]},"black",1);
	print SVG &svg_txt_litc($PaperWidth*1/8+15,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight+$XYTileHeight+$SpaceMarkandXTile+$MarkYStep*($i+1)+10,$MarkFontSize,"black",$LegendFont,$marks[$i]);
	print SVG &svg_rect_litc($PaperWidth*1/2,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight+$XYTileHeight+$SpaceMarkandXTile+$MarkYStep*($i+1),10,10,$colors{$marks[$i+1]},"black",1);
	print SVG &svg_txt_litc($PaperWidth*1/2+15,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight+$XYTileHeight+$SpaceMarkandXTile+$MarkYStep*($i+1)+10,$MarkFontSize,"black",$LegendFont,$marks[$i+1]);

}

print SVG &svg_txt_litc($LeftMargin+$BoxWidth*1/2,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight+$XYTileHeight,$XYTileFontSize,"black",$LegendFont,$Xtile);
print SVG &svg_txt_litc($LeftMargin+$offset*2,$TopMargin+$TileHeight+$TileDownSpace+$BoxHeight/4*2,$XYTileFontSize,"black",$LegendFont,$Ytile,3);

print SVG &svg_end();
close (SVG) ;


# ------------------------------------------------------------------
# Convert to PNG format
# ------------------------------------------------------------------
my $basename = basename($fOut);
my $dirname = dirname($fOut);
chomp(my $cur = `pwd`);
chdir $dirname;
system "$ConvertSVGtoPNG $basename -t pdf";
if (defined $PNG) {
	system("$ConvertSVGtoPNG $basename -w $PaperWidth -h $PaperHeight");
}
chdir $cur;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub svg_paper (){#&svg_paper(width,height,[color])
	my $svg_drawer = "fred_routine"."@"."163\.com";
	chomp $svg_drawer;
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\">\n";
	$line.="<Drawer>$svg_drawer</Drawer>\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub svg_end (){#
	return "</svg>\n";
}

sub svg_txt_litc (){#&svg_txt(x,y,size,color,font,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[6]) {
		$svg_x[6]=0;
	}
	my $svg_matrix='';
	if ($svg_x[6]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[6]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[6]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[6]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"$svg_x[4]\" font-size=\"$svg_x[2]\">$svg_x[5]</text>\n";
	return $line;
}

sub svg_line (){#&svg_line(x1,y1,x2,y2,color,width,[opacity])
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[6]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" opacity=\"$svg_x[6]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_polyline (){#colorfill,colorstroke,width,\@point
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_width=shift(@svg_x);
	my $svg_points=join(" ",@{$svg_x[-1]});
	my $line="<polyline fill=\"$svg_color\" stroke=\"$svg_color2\" stroke-width=\"$svg_width\" points=\"$svg_points\"/>\n";

	return $line;

}

sub svg_rect_litc () {#&svg_rect(x,y,width,height,fillcolor,strokecolor,stroke-width,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[7]) {
		$svg_x[7]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6]\" opacity=\"$svg_x[7]\"/>\n";
	return $line;
}
sub svg_rect_nofill () {#&svg_rect(x,y,width,height,color,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[5]\"/>\n";
	return $line;
}

sub svg_polygon () {#colorfill,colorstroke,coloropacity,point1,point2,...
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_points=join(" ",@svg_x);
	my $line="<polygon fill=\"$svg_color\" stroke=\"$svg_color2\" opacity=\"$svg_trans\" points=\"$svg_points\"/>\n";
	return $line;
}

sub svg_ellipse () {#&svg_ellipse(cx,cy,rx,ry,colorfill,colorstroke,width,[coloropacity])
	my @svg_x=@_;
	my $line= "<ellipse cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" rx=\"$svg_x[2]\" ry=\"$svg_x[3]\" fill=\"$svg_x[4]\" stroke=\"$svg_x[5]\" stroke-width=\"$svg_x[6]\"/>\n";
	if (defined $svg_x[7]) {
		$line="<ellipse cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" rx=\"$svg_x[2]\" ry=\"$svg_x[3]\" fill=\"$svg_x[4]\" stroke=\"$svg_x[5]\" stroke-width=\"$svg_x[6]\" opacity=\"$svg_x[7]\"/>\n";
	}
	return $line;
}

sub svg_circle () {#&svg_circle(cx,cy,r,color)
	my @svg_x=@_;
	my $line="<circle style=\"fill:url(#Gene)\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" r=\"$svg_x[2]\" stroke=\"$svg_x[3]\" stroke-width=\"0\" fill=\"$svg_x[3]\"/>";
	return $line;
}

sub svg_path () {#colorfill,colorstroke,strokewidth,coloropacity,$path
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $width=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_path=shift(@svg_x);
	my $line="<path d= \"$svg_path\" fill=\"$svg_color\" stroke=\"$svg_color2\" stroke-width=\"$width\" opacity=\"$svg_trans\"/>\n";
	return $line;
}

sub LoadData {#
	my ($fIn,$Data)=@_;
	open (IN,"<",$fIn) or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		#my ($class)=(split /\t/,$_)[7];
		my ($class)=$1 if(/.*\s+(\[[A-Z]+\])\s+.*/);
		next if(!defined $class);
		$class=~s/[\[\]]//g;
		foreach my $elem (split "",$class) {
			$$Data{$elem}++;
		}
	}
	close (IN) ;
}
sub OutHash
{
	my ($hash)=@_;
	print OUT "#Category\tID\tClass_Name\tNumbers\n";
	foreach my $cate (sort keys %Category) {
		foreach my $sty (@{$Category{$cate}})
		{
			print OUT "$cate\t"."$sty\t"."$Class{$sty}\t".(${$hash}{$sty}||0)."\n";
		}
	}
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: CogFunClassDrawer.pl
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>
Revised: Yangsh <2011-07-20>
Description:
	Make COG Annotation Classify Graph

Usage:
  Options:
  -i     <file>  input cog annotation result file, forced
  -o     <file>  SVG output file.
  -db    <str>   database of Graph, (COG or KOG)   forced
  -note  <str>   graph head note
  -png       Convert to png format

USAGE
	print $usage;
	exit;
}
