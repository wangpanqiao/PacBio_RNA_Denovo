#!/usr/bin/perl -w
use strict;
use FindBin qw ($Bin);
use File::Basename;
use Getopt::Long;

my $usage =
"
Usage:

  Options:
  -matrix             <file>   rpkm_matrix file.
  -dir                <str>    The project dir.
  -sample             <str>    The prefix of out, (Species Name. etc).
  -fdr                <int>    The FDR threshold value for filter diff-exp result, default is '0.05'
  														  diff-exp filter terms: FDR<0.05 & |FoldChange|>1
  -foldchange         <int>    The FoldChange threshold value for filter diff-exp result, default is '1'
                                diff-exp filter terms: FDR<0.05 & |FoldChange|>1
  -h          Help

For example:
	perl $0 -matrix ../yeast.genes.TMM.fpkm.matrix -dir /RNA_test/denovo_test2/analysis -sample yeast
";

my ($rpkm_matrix,$project_dir,$sample,$fdr,$foldchange,$help);
GetOptions(
	"matrix=s" => \$rpkm_matrix,
	"dir=s" => \$project_dir,
	"sample=s" => \$sample,
	"fdr=i" => \$fdr,
	"foldchange=i" =>\$foldchange,
	"h|help"=>\$help
);
$fdr ||= 0.05;
$foldchange ||= 1;

if (!$rpkm_matrix || !$project_dir || !$sample || $help){
	die "$usage\n";
}
#my $rpkm_matrix = "../yeast.genes.TMM.fpkm.matrix";
#my $project_dir = "/home/leiyoubing/work/work/RNA_test/denovo_test2/analysis";
#my $sample = "yeast";

my %rpkm = ();
open RPKM,"$rpkm_matrix" or die "$!";
my $head = <RPKM>;
chomp($head);
my @head = split /\t/,$head;
my $num = @head;
while(my $line = <RPKM>){
	chomp($line);
	if ($line =~ /^#/){
		next;
	}
	my @arr = split /\t/,$line;
	for (my $i=1;$i<$num;$i++){
		$rpkm{$arr[0]}{$head[$i]} = $arr[$i];
	}
}
close RPKM;

foreach my $diff_result(glob "$project_dir/differential_expression/edgeR.*.dir/*.edgeR.DE_results"){
	my $name = basename $diff_result;
	my $dirname = dirname $diff_result;
	my ($s1,$s2,$prefix);
	if ($name =~ /^(\S+?)_vs_(\S+?).edgeR.DE_results/){
		$s1 = $1;
		$s2 = $2;
	}
	if ($name =~ /^(\S+?_vs_\S+?).edgeR.DE_results/){
		$prefix = $1;
	}
	
	open IN,"$diff_result" or die "$!";
	open OUT,">$dirname/$name.xls" or die "$!";
	open DIFF,">$dirname/$name\_diff_exp.xls" or die "$!";
	<IN>;
	print OUT "#GeneID\tlogFC\tlogCPM\tPValue\tFDR\t$s1\_rpkm\t$s2\_FPKM\n";
	print DIFF "#GeneID\tlogFC\tlogCPM\tPValue\tFDR\t$s1\_rpkm\t$s2\_FPKM\n";
	while(my $line = <IN>){
		chomp($line);
		next if ($line =~ /^\s*$/);
		my @arr = split /\t/,$line;
		print OUT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$rpkm{$arr[0]}{$s1}\t$rpkm{$arr[0]}{$s2}\n";
		my $temp = abs $arr[1];
		if ($arr[4]<$fdr && $temp > $foldchange){
			print DIFF "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$rpkm{$arr[0]}{$s1}\t$rpkm{$arr[0]}{$s2}\n";
		}
	}
	close IN;
	close OUT;
	close DIFF;
	
	unless (-e "$dirname/GO") {mkdir "$dirname/GO", 0755 or die "Cannot make GO: $!";}
	system("Rscript $Bin/topGO_simple.R $dirname/$name\_diff_exp.xls $project_dir/annotation/Result/$sample.Unigenes.fa.GO.list.xls $prefix $dirname/GO");
	
	unless (-e "$dirname/KEGG") {mkdir "$dirname/KEGG", 0755 or die "Cannot make KEGG: $!";}
	unless (-e "$dirname/KEGG/$prefix") {mkdir "$dirname/KEGG/$prefix", 0755 or die "Cannot make KEGG/$prefix: $!";}
	system("perl $Bin/KEGG/kobas_blast.pl -type 1 -gene $dirname/$name\_diff_exp.xls -k $project_dir/annotation/Result/$sample.Unigenes.fa.KEGG.blast.xls -o $dirname/KEGG/$prefix");
}


