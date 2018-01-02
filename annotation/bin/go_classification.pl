#!/usr/bin/perl
#modified by yangying, updata the gene_ontology.1_2.obo to gene_ontology.1_2.0820.obo
use strict;
use warnings;
use lib qw(/home/zhaolili/go-perl-0.15);
use GO::Parser;
die "Usage: Classify genes in the term of gene ontology function.\nExample: perl $0 <gene2go> <output-dir>.\n" unless @ARGV==2;

my $goann=shift;
my $out=shift;
unless(-d $out){
	!system("mkdir $out") or die $!;
}
my $obo="/share/public/database/GO/go_1_2/gene_ontology.1_2.obo";       # obo file
my $parser=new GO::Parser({handler=>'obj',use_cache=>1});
$parser->parse($obo);
my $graph=$parser->handler->graph;

my %exists_go=();
my $it=$graph->create_iterator;
while(my $ni=$it->next_node()){
        $exists_go{$ni->acc}=1;
}



my @root_go=qw(GO:0008150 GO:0005575 GO:0003674);       # root GO terms
my %go_lev1 = ( 'GO:0008150' => 'Biological Process',
                'GO:0005575' => 'Cellular Component',
                'GO:0003674' => 'Molecular Function'
);
my %Namespaces = ( 'biological_process' => 'Biological Process',
                   'cellular_component' => 'Cellular Component',
                   'molecular_function' => 'Molecular Function'
);

my %go_lev2=();
my %go_lev3=();
my %go_lev4=();

my %path=();

foreach my $termLev1(@root_go){
	foreach my $termLev2(@{$graph->get_child_terms($termLev1)}){
		
		$go_lev2{$termLev2->acc}=$termLev2->name;
		foreach my $termLev3(@{$graph->get_child_terms($termLev2)}){
			$go_lev3{$termLev3->acc}=$termLev3->name;
			foreach my $termLev4(@{$graph->get_child_terms($termLev3)}){
				$go_lev4{$termLev4->acc}=$termLev4->name;
				$path{$termLev1}{$termLev2->acc}{$termLev3->acc}{$termLev4->acc}=1;
			}
		}
	}
}

my %gene2go=();
open GO,$goann;
my %genelist=();
while(<GO>){
	chomp;
	my ($gene,@GOs)=split/\t/;
	foreach my $go(@GOs){
		next unless exists($exists_go{$go});
		foreach my $tmp(@{$graph->get_recursive_parent_terms($go)}){
			$gene2go{$tmp->acc}{$gene}=1;
		}
	}
	$genelist{$gene}=1 if @GOs>0;
}
close GO;
my $total_gene=keys %genelist;

open LV2,">$out/GO_classification_count.txt";
print LV2 "## Total annotated genes: $total_gene\n";
print LV2 "#GO ID (Lev2)\tGO Term (Lev2)\tGO Term (Lev1)\tGene Number\n";
#open LV3,">lev3.classification.txt";
open LV4,">$out/GO_classification.xls";
print LV4 "GO ID (Lev1)\tGO Term (Lev1)\tGO ID (Lev2)\tGO Term (Lev2)\tGO ID (Lev3)\tGO Term (Lev3)\tGO ID (Lev4)\tGO Term (Lev4)\tGene Number\tGene List\n";
foreach my $lv1(keys %path){
	foreach my $lv2(keys %{ $path{$lv1}}){
		if(exists $gene2go{$lv2}){
			my $n2=keys %{$gene2go{$lv2}};
	#		my $l2=join(",",keys %{$gene2go{$lv2}});
			print LV2 $lv2."\t".$go_lev2{$lv2}."\t",$go_lev1{$lv1}."\t".$n2."\n";
		}
		foreach my $lv3(keys %{ $path{$lv1}{$lv2}}){
#			if(exists $gene2go{$lv3}){
#				my $n3=keys %{$gene2go{$lv3}};
#				my $l3=join(",",keys %{$gene2go{$lv3}});
#				print LV3 $lv1."\t".$go_lev1{$lv1}."\t".$lv2."\t".$go_lev2{$lv2}."\t".$lv3."\t".$go_lev3{$lv3}."\t".$n3."\t".$l3."\n";
#			}
			foreach my $lv4(keys %{ $path{$lv1}{$lv2}{$lv3}}){
				if(exists $gene2go{$lv4}){
					my $n4=keys %{$gene2go{$lv4}};
					my $l4=join(",",keys %{$gene2go{$lv4}});
					print LV4 $lv1."\t".$go_lev1{$lv1}."\t".$lv2."\t".$go_lev2{$lv2}."\t".$lv3."\t".$go_lev3{$lv3}."\t".$lv4."\t".$go_lev4{$lv4}."\t".$n4."\t".$l4."\n";
				}
			}
		}
	}

}
close LV2;
#close LV3;
close LV4;
