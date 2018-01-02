#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;

my ($info, $seq, $cut, $search, $out,$species,$help);

GetOptions(
	"info=s"	=> \$info,
	"seq=s"		=> \$seq,
	"cut:i"		=> \$cut,
	"search:i"	=> \$search,
	"out=s"		=> \$out,
	"species=s"	=> \$species,
	"help|h"	=> \$help
);

if(!defined($info) || !defined($seq) || !defined($out) || defined($help)){
	exit 0;
}

$cut ||=10;
$search ||=30;
my @suffix = qw(.fa .fasta);
$species ||= basename($seq,@suffix);


my %info=();
$/="\n\/\/\n";
open INFO,$info;
while(<INFO>){
	chomp;
	my ($head,$blast,@array) = split/\n/;
	my @h = split /\t/,$head;
	my @b = split /\t/,$blast;
	my $align_start = $b[1];
	my $align_end = $b[2];

	foreach my $seg(@array){
		my @s = split /\t/,$seg;
		if($s[1] <= $align_start && $s[2] >= $align_end){
			push @{$info{$h[0]}{$h[1]}},$align_start;	##  alignment start
			push @{$info{$h[0]}{$h[1]}},$align_end;		##  alignment end
			push @{$info{$h[0]}{$h[1]}},$s[1];		##  proper segment start
			push @{$info{$h[0]}{$h[1]}},$s[2];		##  proper segment end
		}
	}
}
close INFO;

unless(-d $out){
	!system("mkdir $out") or die $!;
}

my $fa_obj=Bio::SeqIO->new(-file=>$seq,-format=>"Fasta");
my $cds_obj=Bio::SeqIO->new(-file=>">$out/$species.cds.fasta",-format=>'Fasta');
my $pep_obj=Bio::SeqIO->new(-file=>">$out/$species.pep.fasta",-format=>'Fasta');
my $utr5_obj=Bio::SeqIO->new(-file=>">$out/$species.5_UTR.fasta",-format=>'Fasta');
my $utr3_obj=Bio::SeqIO->new(-file=>">$out/$species.3_UTR.fasta",-format=>'Fasta');

while(my $gene = $fa_obj->next_seq()){
	my $id=$gene->display_id();
	if(exists($info{$id})){
		my $count=1;
		my ($cds_start,$cds_end);
		my $hln = keys %{ $info{$id}};
		if($hln>1){
			next;
		}
		foreach my $f(keys %{ $info{$id}}){
			my ($as,$ae,$ss,$se) = @{$info{$id}{$f}};
			my $gln = $gene->length;
			#########################################################################
			my $show_id = $id;
			$id.="_orf$count" if $count>1;
			my $desc = "\tstart:".$as."\t"."end:".$se."\tframe:$f";
			my $utr5_desc = "\tstart:1"."\t"."end:".($as-1);
			my $utr3_desc = "\tstart:".($se+1)."\t"."end:".$gln;
			#########################################################################
			if($f<0){
				$gene = $gene->revcom;

				$desc = "\tstart:".$ae."\t"."end:".$ss."\tframe:$f";
				$utr5_desc = "\tstart:".$gln."\t"."end:".($ae+1);
				$utr3_desc = "\tstart:".$ss."\t"."end:1";

				my $tmp = $gln-$ae+1;
				$ae = $gln-$as+1;
				$as = $tmp;
				$tmp = $gln-$se+1;
				$se = $gln-$ss+1;
				$ss = $tmp;
			}
			$cds_start = $as;
			$cds_end = $se;

			my ($cds, $pep, $utr5, $utr3);
			if(abs($as-$ss)<3){
				$cds = $gene->trunc($cds_start,$cds_end);
				$pep = $cds->translate();
				$cds->display_id($show_id);
				$cds->desc($desc."\tlength:".$cds->length());
				$pep->display_id($show_id);
				$pep->desc($desc."\tlength:".$pep->length());
				$cds_obj->write_seq($cds);
				$pep_obj->write_seq($pep);
				
				if($gln-$cds_end>$cut){
					$utr3 = $gene->trunc($cds_end+1,$gln);
					$utr3->display_id($show_id);
					$utr3->description($utr3_desc."\tlength:".$utr3->length());
					$utr3_obj->write_seq($utr3);
				}
				next;
			}
			###########################################
			{
				my $ts = $cds_start-$search;
				my $te = $cds_start+$search;
				$te = $gln if $te>$gln;
				$ts=abs($f) if $ts<=0;
				my $motif = $gene->subseq($ts,$te);
				my $atg=&searchATG($motif);
				unless($atg == 100000){
					$cds_start = $ts+$atg;
					$utr5_desc = "\tstart:1"."\t"."end:".($cds_start+1);
					$utr5_desc = "\tstart:$gln"."\tend:".($gln-$cds_start+1) if $f<0;
				}
			}
			###########################################
			$count++;
			$cds = $gene->trunc($cds_start,$cds_end);
			$pep = $cds->translate();
			$cds->display_id($show_id);
			$cds->desc($desc."\tlength:".$cds->length());
			$pep->display_id($show_id);
			$pep->desc($desc."\tlength:".$pep->length());
			$cds_obj->write_seq($cds);
			$pep_obj->write_seq($pep);
			
			if($cds_start>$cut){
				$utr5 = $gene->trunc(1,$cds_start-1);
				$utr5->display_id($show_id);
				$utr5->description($utr5_desc."\tlength:".$utr5->length());
				$utr5_obj->write_seq($utr5);
			}
			if($gln-$cds_end>$cut){
				$utr3 = $gene->trunc($cds_end+1,$gln);
				$utr3->display_id($show_id);
				$utr3->description($utr3_desc."\tlength:".$utr3->length());
				$utr3_obj->write_seq($utr3);
			}
		}
	}
}

my %stop=("TAA"=>1,"TGA"=>1,"TAG"=>1);

sub searchATG{
	my $str =shift @_;
	my $len = length($str)/3;
	my $re = 100000;
	my $first=1;
	foreach my $i(1..$len){
		my $s = $i*3-1;
		my $e = $i*3;
		my $tmp = substr($str, $s, $e);
		if($tmp eq "ATG" && $first){
			$re = $s;
			$first=0;
		#	last;
		}
		elsif(exists($stop{$tmp})){
			$re = 100000;
			$first=1;
		}
	}
	return $re;
}

