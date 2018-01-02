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
  -repeat             <file>   The repeat condition-file.
	-fdr                <int>    The FDR threshold value for filter diff-exp result, default is '0.05'
  														  diff-exp filter terms: FDR<0.05 & |FoldChange|>1
  -foldchange         <int>    The FoldChange threshold value for filter diff-exp result, default is '1'
                                diff-exp filter terms: FDR<0.05 & |FoldChange|>1
  -h          Help

For example:
	perl $0 -matrix ../yeast.genes.TMM.fpkm.matrix -dir /RNA_test/denovo_test2/analysis -sample yeast -repeat samples_condition_file
";

my ($rpkm_matrix,$project_dir,$sample,$repeat,$fdr,$foldchange,$help);
GetOptions(
	"matrix=s" => \$rpkm_matrix,
	"dir=s" => \$project_dir,
	"sample=s" => \$sample,
	"repeat=s" => \$repeat,
	"fdr=f" => \$fdr,
	"foldchange=i" =>\$foldchange,
	"h|help"=>\$help
);
$fdr ||= 0.05;
$foldchange ||= 1;

if (!$rpkm_matrix || !$project_dir || !$sample || !$repeat || $help){
	die "$usage\n";
}
#my $rpkm_matrix = "../yeast.genes.TMM.fpkm.matrix";
#my $project_dir = "/home/leiyoubing/work/work/RNA_test/denovo_test2/analysis";
#my $sample = "yeast";

my (%condition,%cond_tmp) = ();
open REPEAT,"$repeat" || die "$!";
while(my $line = <REPEAT>){
	chomp($line);
	my @tmp = split /\t/,$line;
	$condition{$tmp[1]} = $tmp[0];
	$cond_tmp{$tmp[0]} .= "$tmp[1],";
}
close REPEAT;

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

foreach my $diff_result(glob "$project_dir/edgeR.*.dir/*.edgeR.DE_results"){
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
	
	my $condA = $cond_tmp{$s1};
	my $condB = $cond_tmp{$s2};
	$condA =~ s/,$//;
	$condB =~ s/,$//;
	my @A = split /,/,$condA;
	my @B = split /,/,$condB;
	
	open IN,"$diff_result" or die "$!";
	open OUT,">$dirname/$name.xls" or die "$!";
	open DIFF,">$dirname/$name\_diff_exp.xls" or die "$!";
	<IN>;
	print OUT "#GeneID\tlogFC\tlogCPM\tPValue\tFDR";
	print DIFF "#GeneID\tlogFC\tlogCPM\tPValue\tFDR";
	for(my $i=0;$i<@A;$i++){
		print OUT "\t$A[$i]\_FPKM_$condition{$A[$i]}";
		print DIFF "\t$A[$i]\_FPKM_$condition{$A[$i]}";
	}
	for(my $j=0;$j<@B;$j++){
		print OUT "\t$B[$j]\_FPKM_$condition{$B[$j]}";
		print DIFF "\t$B[$j]\_FPKM_$condition{$B[$j]}";
	}
	print OUT "\n";
	print DIFF "\n";
	
	while(my $line = <IN>){
		chomp($line);
		next if ($line =~ /^\s*$/);
		my @arr = split /\t/,$line;
		print OUT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]";
		for(my $i=0;$i<@A;$i++){
			print OUT "\t$rpkm{$arr[0]}{$A[$i]}";
		}
		for(my $j=0;$j<@B;$j++){
			print OUT "\t$rpkm{$arr[0]}{$B[$j]}";
		}
		print OUT "\n";
		
		my $temp = abs $arr[1];
		if ($arr[4]<$fdr && $temp > $foldchange){
			print DIFF "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]";
			for(my $i=0;$i<@A;$i++){
				print DIFF "\t$rpkm{$arr[0]}{$A[$i]}";
			}
			for(my $j=0;$j<@B;$j++){
				print DIFF "\t$rpkm{$arr[0]}{$B[$j]}";
			}
			print DIFF "\n";
		}
	}
	close IN;
	close OUT;
	close DIFF;

	unless (-e "$dirname/GO") {mkdir "$dirname/GO", 0755 or die "Cannot make GO: $!";}
	system("Rscript $Bin/topGO_simple.R $dirname/$name\_diff_exp.xls $project_dir/RefIso.fa.GO.list.xls $prefix $dirname/GO");
	
	unless (-e "$dirname/KEGG") {mkdir "$dirname/KEGG", 0755 or die "Cannot make KEGG: $!";}
	unless (-e "$dirname/KEGG/$prefix") {mkdir "$dirname/KEGG/$prefix", 0755 or die "Cannot make KEGG/$prefix: $!";}
	system("perl $Bin/KEGG/kobas_blast.pl -type 1 -gene $dirname/$name\_diff_exp.xls -k $project_dir/RefIso.fa.KEGG.blast.xls -o $dirname/KEGG/$prefix");
}
