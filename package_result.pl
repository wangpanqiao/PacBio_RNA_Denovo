#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin '$Bin';
use File::Find;
use Getopt::Long;
use Cwd qw (abs_path);

if (@ARGV < 3){
	die "perl $0 analysis_path prefix_of_sample pacbio_sample_list\n";
}
my ($path,$prefix,$paclist) = @ARGV;
$paclist=abs_path($paclist);
#####################################################################################
#             make home directionary for the web report result                      #
#####################################################################################
unless (-e "$path/PacBio_RNA_denovo_report") {mkdir "$path/PacBio_RNA_denovo_report", 0755;}
unless (-e "$path/PacBio_RNA_denovo_report/result") {mkdir "$path/PacBio_RNA_denovo_report/result", 0755;}
unless (-e "$path/PacBio_RNA_denovo_report/result/primary") {mkdir "$path/PacBio_RNA_denovo_report/result/primary", 0755;}
if (-e "$path/result/IsoSeq" && !-e "$path/PacBio_RNA_denovo_report/result/primary/IsoSeq") {mkdir "$path/PacBio_RNA_denovo_report/result/primary/IsoSeq", 0755;}
if (-e "$path/result/CD-hit-est" && !-e "$path/PacBio_RNA_denovo_report/result/primary/CD-hit-est") {mkdir "$path/PacBio_RNA_denovo_report/result/primary/CD-hit-est", 0755;}
if (-e "$path/result/Transdecoder" && !-e "$path/PacBio_RNA_denovo_report/result/primary/Transdecoder") {mkdir "$path/PacBio_RNA_denovo_report/result/primary/Transdecoder", 0755;}
if (-e "$path/result/OrthoMCL" && !-e "$path/PacBio_RNA_denovo_report/result/primary/OrthoMCL") {mkdir "$path/PacBio_RNA_denovo_report/result/primary/OrthoMCL", 0755;}
if (-e "$path/result/Annotation" && !-e "$path/PacBio_RNA_denovo_report/result/primary/Annotation") {mkdir "$path/PacBio_RNA_denovo_report/result/primary/Annotation", 0755;}
if (-e "$path/result/SSR" && !-e "$path/PacBio_RNA_denovo_report/result/primary/SSR") {mkdir "$path/PacBio_RNA_denovo_report/result/primary/SSR", 0755;}
if (-e "$path/result/LncRNA" && !-e "$path/PacBio_RNA_denovo_report/result/primary/LncRNA") {mkdir "$path/PacBio_RNA_denovo_report/result/primary/LncRNA", 0755;}
unless (-e "$path/PacBio_RNA_denovo_report/result/advanced"){
if (-e "$path/result/LSC" || -e "$path/result/RSEM" || -e "$path/result/edgeR") {mkdir "$path/PacBio_RNA_denovo_report/result/advanced", 0755;}}
if (-e "$path/result/LSC" && !-e "$path/PacBio_RNA_denovo_report/result/advanced/LSC") {mkdir "$path/PacBio_RNA_denovo_report/result/advanced/LSC", 0755;}
if (-e "$path/result/RefIso" && !-e "$path/PacBio_RNA_denovo_report/result/advanced/RefIso") {mkdir "$path/PacBio_RNA_denovo_report/result/advanced/RefIso", 0755;}
if (-e "$path/result/RSEM" && !-e "$path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation") {mkdir "$path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation", 0755;}
if (-e "$path/result/edgeR" && ! -e "$path/PacBio_RNA_denovo_report/result/advanced/Differential_expression") {mkdir "$path/PacBio_RNA_denovo_report/result/advanced/Differential_expression", 0755;}
system("cp $path/result/base_info.json $path/PacBio_RNA_denovo_report/");
system("cp $path/PacBio_sample.xls $path/PacBio_RNA_denovo_report/result/primary/IsoSeq");
#####################################################################################
#                              collect the output of IsoSeq                         #
#####################################################################################
open (IN,"$paclist");
while (my $line=<IN>){
chomp $line;
my @array=split /\t/,$line;
my $sample='';
my $cDNA='';
my @cDNAs=();
my @frag=("under1k","between1k2k","between2k3k","above3k","between1k4k","between3k6k","between5k10k");
my $dir1="$path/result/IsoSeq";
my $dir2="$path/PacBio_RNA_denovo_report/result/primary/IsoSeq";
foreach my $ele (@frag){
if (-e "$path/result/IsoSeq/$array[0]/$array[0]\_$ele"){
push @cDNAs,$ele;
unless (-e "$dir2/$array[0]") {mkdir "$dir2/$array[0]", 0755;}
unless (-e "$dir2/$array[0]/$array[0]\_$ele") {mkdir "$dir2/$array[0]/$array[0]\_$ele", 0755;}
unless (-e "$dir2/$array[0]/$array[0]\_$ele/results") {mkdir "$dir2/$array[0]/$array[0]\_$ele/results", 0755;}
system("cp $dir1/$array[0]/$array[0]\_$ele/PacBio_rawdata*.png $dir2/$array[0]/$array[0]\_$ele/results");
system("cp $dir1/$array[0]/$array[0]\_$ele/results/*.json $dir2/$array[0]/$array[0]\_$ele/results");
system("cp $dir1/$array[0]/$array[0]\_$ele/results/roi_readlength_hist.png $dir2/$array[0]/$array[0]\_$ele/results/$array[0]\_$ele\_roi_readlength_hist.png");
system("cp $dir1/$array[0]/$array[0]\_$ele/results/roi_accuracy_hist.png $dir2/$array[0]/$array[0]\_$ele/results/$array[0]\_$ele\_roi_accuracy_hist.png");
system("cp $dir1/$array[0]/$array[0]\_$ele/results/roi_npasses_hist.png $dir2/$array[0]/$array[0]\_$ele/results/$array[0]\_$ele\_roi_npasses_hist.png");
system("cp $dir1/$array[0]/$array[0]\_$ele/results/classify_summary.png $dir2/$array[0]/$array[0]\_$ele/results/$array[0]\_$ele\_classify_summary.png");
system("cp $dir1/$array[0]/$array[0]\_$ele/results/classify_summary.pdf $dir2/$array[0]/$array[0]\_$ele/results/$array[0]\_$ele\_classify_summary.pdf");
system("cp $dir1/$array[0]/$array[0]\_$ele/results/*length.xls $dir2/$array[0]/$array[0]\_$ele/results");
system("cp $dir1/$array[0]/$array[0]\_$ele/data/reads_of_insert.fast* $dir2/$array[0]/$array[0]\_$ele/results");
system("cp $dir1/$array[0]/$array[0]\_$ele/data/isoseq*fl*fasta $dir2/$array[0]/$array[0]\_$ele/results");
system("cp $dir1/$array[0]/$array[0]\_$ele/data/*consensus_isoforms* $dir2/$array[0]/$array[0]\_$ele/results");
}
}
$cDNA=join (",",@cDNAs);
#system("python /share/work2/liuying/pipeline/PacBio_RNA_Denovo/Bin/IsoSeq/iso_stat_sum.py --indir $dir2/ --sample $array[0] --cDNA $cDNA --outdir $dir2/$array[0]");


system("cat $dir2/$array[0]/*/results/roi_length.xls >$dir2/$array[0]/$array[0]\_roi_length.xls");
system("Rscript $Bin/IsoSeq/plot_length.R $dir2/$array[0]/$array[0]\_roi_length.xls $dir2/$array[0]/$array[0]\_roi_length_density");
system("cat $dir2/$array[0]/*/results/flnc_length.xls >$dir2/$array[0]/$array[0]\_flnc_length.xls");
system("Rscript $Bin/IsoSeq/plot_length.R $dir2/$array[0]/$array[0]\_flnc_length.xls $dir2/$array[0]/$array[0]\_flnc_length_density");
system("cat $dir2/$array[0]/*/results/consensus_length.xls >$dir2/$array[0]/$array[0]\_consensus_length.xls");
system("Rscript $Bin/IsoSeq/plot_length.R $dir2/$array[0]/$array[0]\_consensus_length.xls $dir2/$array[0]/$array[0]\_consensus_length_density");
system("rm -r $dir2/$array[0]/*/results/*.json");
system("rm -r $dir2/$array[0]/*/results/*length.xls");
system("rm -r $dir2/$array[0]/*length.xls");
system("python $Bin/IsoSeq/iso_stat_sum.py --indir $dir1/ --sample $array[0] --cDNA $cDNA --outdir $dir2/$array[0]");
system("$Bin/IsoSeq/merge_roi_mulsamples_data --indir $dir2/ -name roi_stat.xls --sample $array[0] --outfile $dir2/all_roi_stat.xls");
system("$Bin/IsoSeq/merge_roi_mulsamples_data --indir $dir2/ -name classify_stat.xls --sample $array[0] --outfile $dir2/all_classify_stat.xls");
system("$Bin/IsoSeq/merge_roi_mulsamples_data --indir $dir2/ -name cluster_stat.xls --sample $array[0] --outfile $dir2/all_cluster_stat.xls");


#####################################################################################
#                              collect the output of LSC                            #
#####################################################################################
if (-e "$path/result/LSC/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/advanced/LSC/$array[0]");
system("cp  $path/result/LSC/$array[0]/cat_lq_consensus_isoforms.fasta $path/PacBio_RNA_denovo_report/result/advanced/LSC/$array[0]/pre_LSC.fasta");
system("cp  $path/result/LSC/$array[0]/full_LR.fa $path/PacBio_RNA_denovo_report/result/advanced/LSC/$array[0]/LSC.fasta");
}
#####################################################################################
#                       collect the output of CD-hit-est                            #
#####################################################################################
if (-e "$path/result/CD-hit-est/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/CD-hit-est/$array[0]");
system("cp  $path/result/CD-hit-est/$array[0]/cat_consensus_isoforms.fasta $path/PacBio_RNA_denovo_report/result/primary/CD-hit-est/$array[0]/pre_CD-hit-est.fasta");
system("cp  $path/result/CD-hit-est/$array[0]/UniIso.fa $path/PacBio_RNA_denovo_report/result/primary/CD-hit-est/$array[0]/CD-hit-est.fasta");
system("cp  $path/result/CD-hit-est/$array[0]/*.fa.clstr $path/PacBio_RNA_denovo_report/result/primary/CD-hit-est/$array[0]/CD-hit-est.clstr.txt");
}
#####################################################################################
#                       collect the output of Transdecoder                          #
#####################################################################################
if (-e "$path/result/Transdecoder/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/Transdecoder/$array[0]");
system("cp  $path/result/Transdecoder/$array[0]/UniIso.fa.transdecoder.gff3 $path/PacBio_RNA_denovo_report/result/primary/Transdecoder/$array[0]/transdecoder.gff3");
system("cp  $path/result/Transdecoder/$array[0]/UniIso.fa.transdecoder.cds.txt $path/PacBio_RNA_denovo_report/result/primary/Transdecoder/$array[0]/transdecoder.cds.txt");
system("cp  $path/result/Transdecoder/$array[0]/UniIso.fa.transdecoder.pep.txt $path/PacBio_RNA_denovo_report/result/primary/Transdecoder/$array[0]/transdecoder.pep.txt");
}
#####################################################################################
#                       collect the output of OrthoMCL                              #
#####################################################################################
if (-e "$path/result/OrthoMCL/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/OrthoMCL/$array[0]");
system("cp  $path/result/OrthoMCL/$array[0]/$prefix.group.xls $path/PacBio_RNA_denovo_report/result/primary/OrthoMCL/$array[0]/group.xls");
}
#####################################################################################
#                       collect the output of Annotation                            #
#####################################################################################
if (-e "$path/result/Annotation/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]");
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/GO");
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KEGG");
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KOG");
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/NT_NR_SW");
system("cp $path/result/Annotation/$array[0]/Result/Function_Annotation.stat.xls $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]");
system("cp $path/result/Annotation/$array[0]/Result/Integrated_Function.annotation.xls $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]");
system("cp $path/result/Annotation/$array[0]/Result/UniIso.fa.Anno_Venn.png $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/Anno_Venn.png");
system("cp $path/result/Annotation/$array[0]/Result/GO_classification* $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/GO");
system("cp $path/result/Annotation/$array[0]/Result/UniIso.fa.GO.anno.xls $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/GO/GO.anno.xls");
system("cp $path/result/Annotation/$array[0]/Result/KEGG_classification_count.txt $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KEGG");
system("sed -i '1d' $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KEGG/KEGG_classification_count.txt");
system("cp $path/result/Annotation/$array[0]/Result/KEGG_classification.pdf $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KEGG");
system("cp $path/result/Annotation/$array[0]/Result/KEGG_classification.png $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KEGG");
system("cp $path/result/Annotation/$array[0]/Result/UniIso.fa.KEGG.gene2ko.xls $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KEGG/KEGG.gene2ko.xls");
system("cp $path/result/Annotation/$array[0]/Result/*KOG*class* $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/KOG");
system("cp $path/result/Annotation/$array[0]/Result/UniIso.fa.NR.anno.xls $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/NT_NR_SW/NR.anno.xls");
system("cp $path/result/Annotation/$array[0]/Result/UniIso.fa.NT.anno.xls $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/NT_NR_SW/NT.anno.xls");
system("cp $path/result/Annotation/$array[0]/Result/UniIso.fa.Swissprot.anno.xls $path/PacBio_RNA_denovo_report/result/primary/Annotation/$array[0]/NT_NR_SW/Swissprot.anno.xls");
}
#####################################################################################
#                       collect the output of SSR                                   #
#####################################################################################
if (-e "$path/result/SSR/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/SSR/$array[0]");
system("cp  $path/result/SSR/$array[0]/UniIso.fa.name.convert $path/PacBio_RNA_denovo_report/result/primary/SSR/$array[0]");
system("cp  $path/result/SSR/$array[0]/UniIso.fa.results.xls $path/PacBio_RNA_denovo_report/result/primary/SSR/$array[0]");
system("cp  $path/result/SSR/$array[0]/ssr_density.* $path/PacBio_RNA_denovo_report/result/primary/SSR/$array[0]");
system("cp  $path/result/SSR/$array[0]/Transcripts.fa.SSR_motifs_distribution.* $path/PacBio_RNA_denovo_report/result/primary/SSR/$array[0]");
}
#####################################################################################
#                       collect the output of LncRNA                                #
#####################################################################################
if (-e "$path/result/LncRNA/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/primary/LncRNA/$array[0]");
system("cp  $path/result/LncRNA/$array[0]/UniIso_ncrna.fa $path/PacBio_RNA_denovo_report/result/primary/LncRNA/$array[0]/ncrna.fa");
system("cp  $path/result/LncRNA/$array[0]/lncRNA_length_cout.png $path/PacBio_RNA_denovo_report/result/primary/LncRNA/$array[0]");
}
}
close IN;
#####################################################################################
#                       collect the output of RefIso                                #
#####################################################################################
if (-e "$path/result/RefIso"){
system("cp $path/result/RefIso/RefIso.fa $path/PacBio_RNA_denovo_report/result/advanced/RefIso");
if (-e "$path/result/RefIso/Multiple_isoforms.fasta"){
system("cp $path/result/RefIso/Multiple_isoforms.fasta $path/PacBio_RNA_denovo_report/result/advanced/RefIso/pre_RefIso.fa");}
}
#####################################################################################
#                       collect the output of expression evaluation                 #
#####################################################################################
if (-e "$path/result/RSEM"){
open SH,"$path/result/RSEM/sample_list.txt";
my $line=<SH>;
my $combine="";
chomp $line;
my @sample=split /,/,$line;
foreach my $ele(@sample){$combine=$combine." ".$path."/result/RSEM/".$ele."/".$ele.".Readcount_FPKM.xls";}
system("$Bin/expression/merge_FPKM_single_table $combine >$path/result/RSEM/merged.fpkm");
system("$Bin/expression/merge_readcounts_single_table $combine >$path/result/RSEM/merged.readcount");
system("$Bin/expression/calrowmeans -rpkm $path/result/RSEM/merged.fpkm -group $line -groupname $line -out-rowmeans $path/result/RSEM/rowmeans_fpkm.xls");
system("$Bin/expression/plot.FPKM_v3.1 -r $path/result/RSEM/rowmeans_fpkm.xls -fpkm $path/result/RSEM/merged.fpkm -rc $path/result/RSEM/merged.readcount -len $path/result/RSEM/$sample[0]/geneINFO -output $path/result/RSEM");
my $dir = "$path/result/RSEM";
my @files=();
sub wanted {
    if ( -d $File::Find::dir ) {
        if ( $File::Find::dir =~ /^$path\/result\/RSEM\/([^\/]+)$/ )#give the standard
 {
            push @files,$File::Find::dir."\t".$1;
        }
        }
 return(@files);
        }
 find (\&wanted,$dir);#find the directory within $dir matching the stardard, and return one array @files. 
my %count = ();
my @uniq = grep {!$count{$_} ++} @files;#uniq the array
foreach my $each (@uniq){
my @newpath=split /\t/,$each;
system("mkdir  $path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation/$newpath[1]");
system("cp  $path/result/RSEM/FPKM_boxplot* $path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation/");
system("cp  $path/result/RSEM/FPKM_density_distribution* $path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation/");
system("cp  $newpath[0]/*TranscriptExp.xls $path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation/$newpath[1]");
system("cp  $newpath[0]/*.mapped.stat.xls $path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation/$newpath[1]");
}
system("cat $path/result/RSEM/*/*.mapped.stat.2.xls | sed '1i Sample	Total Reads	Mapped Reads(ratio)	Uniq mapped Reads(ratio)	Multi mapped Reads(ratio)' > $path/PacBio_RNA_denovo_report/result/advanced/Expression_evaluation/Total_mapped.stat.xls");
}
#####################################################################################
#                       collect the output of differential expression               #
#####################################################################################
if (-e "$path/result/edgeR"){
#	system("cp -r $path/result/edgeR/* $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("cp $path/result/edgeR/edgeR.*.dir/*edgeR.DE_results_diff_exp.xls $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("cp $path/result/edgeR/edgeR.*.dir/*edgeR.DE_results.MA* $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("cp $path/result/edgeR/edgeR.*.dir/*edgeR.DE_results.Volcano* $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("cp $path/result/edgeR/edgeR.*.dir/*genes_vs_samples_heatmap.png $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("cp $path/result/edgeR/edgeR.*.dir/*sample_cor_matrix.png $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("cp -r $path/result/edgeR/edgeR.*.dir/GO $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("cp -r $path/result/edgeR/edgeR.*.dir/KEGG/ $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression");
	system("rm -r $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression/KEGG/*/*.annotate");
	system("rm -r $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression/KEGG/*/*.blast");
	system("rm -r $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression/KEGG/*/*.identify");
	system("rm -r $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression/KEGG/*/*.Kegg.ko2gene.xls");
	system("rm -r $path/PacBio_RNA_denovo_report/result/advanced/Differential_expression/KEGG/*/*.tmp.xls");
}
#####################################################################################
#                                  produce the web report                           #
#####################################################################################
	open SH, ">$path/PacBio_RNA_denovo_report/web.sh";
	print SH "/share/public/software/node-v6.11.1/bin/node $Bin/jade2html/pac_rna_denovo.js ./";
	close SH;
	system("cd $path/PacBio_RNA_denovo_report/ && qsub -cwd web.sh");
