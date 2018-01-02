#!/usr/bin/perl -w
use strict;
use LWP::Simple;
use File::Basename;
use Cwd qw(abs_path realpath getcwd);

my ($kos2gene,$pathway_enrichment,$outdir)=@ARGV;

unless(@ARGV==3){
        print "\n\tUsage:\n\t\tperl $0 <kos2gene file> <pathway enrichment file> <outdir>\n";
        print "\n\tExample:\n\t\tperl $0 Up_kos2gene.xls Up_pathway_enrichment.xls /*/up_pathway_map\n";
	print "\n\tNote: outdir need to designate absolute path !\n\n"; 	
        exit;
}

my %kos2gene;
###load kos2gene 
open KOS2GENE,$kos2gene or die;
<KOS2GENE>;
while(<KOS2GENE>){
	chomp;
	my @arr=split /\t/;
	my ($ko_id,$ko_name)=$arr[0]=~/(K\d+)\((.+)\)/;
	push @{$kos2gene{$ko_id}},($ko_name,$arr[1]);
}
close KOS2GENE;

###load pathway enrichment file
my $file=basename $pathway_enrichment;   
open PATHWAY,$pathway_enrichment or die;
<PATHWAY>;
while(<PATHWAY>){
	chomp;
	my ($ko,$url)=(split /\t/,$_)[2,8];
	my $content=get($url);
	next unless $content ;

	######extract content of map position###############
	my ($lines)=$content=~/(<map name="mapdata">.*?<\/map>)/s;
	my @lines=split /\n/,$lines;
	my $head=html_head($ko);
	my $body=html_body($ko);

	###produce map.html###
	system "mkdir $outdir" unless ( -d $outdir );
	open HTML,">$outdir/$ko.html" or die;
	print HTML "$head\n<body>\n";
	foreach (@lines){
		if ($_ eq "<map name=\"mapdata\">"){
			print HTML "<map name=\"$ko\">\n";
			next;
		}
		
		$_=~s/href="/href="http:\/\/www.genome.jp/;	
		$_=~s/\ttitle=.* \/>/ \/>/;
		###onmouseover###
		my @onmouseover;
		my @ko=$_=~/(K\d+)/g;	
		foreach (@ko){
		 	push @onmouseover,"<li>$_($kos2gene{$_}->[0]): $kos2gene{$_}->[1]</li>" if $kos2gene{$_};
		}

		if(@onmouseover>0){
			my $onmouseover=join "",@onmouseover;
			$onmouseover="onmouseover=\'javascript: showInfo(\"<ul><li style=\\\"color: #f00;\\\">Gene<ul>$onmouseover</ul></li></ul>\");\'";
			$_=~s/ \/>/\t$onmouseover \/>/;
		}
		print HTML "$_\n";
	}
	
	###html body###
	print HTML "$body";
	###download png###
	if($content=~/<img src=\"(.+png)\"/){
    		system("wget -q http://www.genome.jp/$1 -P $outdir");
	}
	
}
close HTML;
close PATHWAY;

######sub function######
sub html_head{
	my $map=shift;
	my $head="<html>
<head>
<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">
<title>$map</title>
<style type=\"text/css\">
<!--
area {cursor: pointer;}
-->
</style>
<script type=\"text/javascript\">
<!--
function showInfo(info) {
        obj = document.getElementById(\"result\");
        obj.innerHTML = \"<div style=\'cursor: pointer; position: absolute; right: 5px; color: #000;\' onclick=\'javascript: document.getElementById(\\\"result\\\").style.display = \\\"none\\\";\' title=\'close\'>X</div>\" + info;
        obj.style.top = document.body.scrollTop;
        obj.style.left = document.body.scrollLeft;
        obj.style.display = \"\";
}
//->
</script>
</head>";
	return($head);
}

sub html_body {
	my $map=shift;
	my $body="<img src=\'./$map.png\' usemap=\'#$map\' />
<div id=\'result\' style=\'position: absolute; width: 50%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;\' onmouseover=\"javascript: this.style.filter = \'alpha(opacity=100)\'; this.style.opacity = 1;\" onmouseout=\"javascript: this.style.filter = \'alpha(opacity=95)\'; this.style.opacity = 0.95;\"></div>
</body></html>";
}

