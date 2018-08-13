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
if (-e "$path/IsoSeq" && !-e "$path/PacBio_RNA_denovo_report/result/IsoSeq") {mkdir "$path/PacBio_RNA_denovo_report/result/IsoSeq", 0755;}
if (-e "$path/CD-hit-est" && !-e "$path/PacBio_RNA_denovo_report/result/CD-hit-est") {mkdir "$path/PacBio_RNA_denovo_report/result/CD-hit-est", 0755;}
if (-e "$path/Transdecoder" && !-e "$path/PacBio_RNA_denovo_report/result/Transdecoder") {mkdir "$path/PacBio_RNA_denovo_report/result/Transdecoder", 0755;}
if (-e "$path/OrthoMCL" && !-e "$path/PacBio_RNA_denovo_report/result/OrthoMCL") {mkdir "$path/PacBio_RNA_denovo_report/result/OrthoMCL", 0755;}
if (-e "$path/Annotation" && !-e "$path/PacBio_RNA_denovo_report/result/Annotation") {mkdir "$path/PacBio_RNA_denovo_report/result/Annotation", 0755;}
if (-e "$path/SSR" && !-e "$path/PacBio_RNA_denovo_report/result/SSR") {mkdir "$path/PacBio_RNA_denovo_report/result/SSR", 0755;}
if (-e "$path/LncRNA" && !-e "$path/PacBio_RNA_denovo_report/result/LncRNA") {mkdir "$path/PacBio_RNA_denovo_report/result/LncRNA", 0755;}
if (-e "$path/LSC" && !-e "$path/PacBio_RNA_denovo_report/result/LSC") {mkdir "$path/PacBio_RNA_denovo_report/result/LSC", 0755;}
if (-e "$path/RefIso" && !-e "$path/PacBio_RNA_denovo_report/result/RefIso") {mkdir "$path/PacBio_RNA_denovo_report/result/RefIso", 0755;}
if (-e "$path/RSEM" && !-e "$path/PacBio_RNA_denovo_report/result/Expression_evaluation") {mkdir "$path/PacBio_RNA_denovo_report/result/Expression_evaluation", 0755;}
if (-e "$path/edgeR" && ! -e "$path/PacBio_RNA_denovo_report/result/Differential_expression") {mkdir "$path/PacBio_RNA_denovo_report/result/Differential_expression", 0755;}
if (-e "$path/AS" && ! -e "$path/PacBio_RNA_denovo_report/result/AS") {mkdir "$path/PacBio_RNA_denovo_report/result/AS", 0755;}
system("cp $path/base_info.json $path/PacBio_RNA_denovo_report/");
system("cp $path/PacBio_sample.xls $path/PacBio_RNA_denovo_report/result/IsoSeq");
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
my @frag=("0-1k","1-2k","2-3k","3-4k","1-4k","1-6k","3-6k","3-10k","4-10k","5-10k");
my $dir1="$path/IsoSeq";
my $dir2="$path/PacBio_RNA_denovo_report/result/IsoSeq";
foreach my $ele (@frag){
if (-e "$path/IsoSeq/$array[0]/$array[0]\_$ele"){
push @cDNAs,$ele;
unless (-e "$dir2/$array[0]") {mkdir "$dir2/$array[0]", 0755;}
unless (-e "$dir2/$array[0]/$array[0]\_$ele") {mkdir "$dir2/$array[0]/$array[0]\_$ele", 0755;}
## system("cp $dir1/$array[0]/$array[0]\_$ele/PacBio_rawdata*.png $dir2/$array[0]/$array[0]\_$ele/results");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/*.json $dir2/$array[0]/$array[0]\_$ele");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/*_readlength_hist.png $dir2/$array[0]/$array[0]\_$ele/$array[0]\_$ele\_roi_readlength_hist.png");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/*_accuracy_hist.png $dir2/$array[0]/$array[0]\_$ele/$array[0]\_$ele\_roi_accuracy_hist.png");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/*_npasses_hist.png $dir2/$array[0]/$array[0]\_$ele/$array[0]\_$ele\_roi_npasses_hist.png");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/*_classify_summary.* $dir2/$array[0]/$array[0]\_$ele");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/*length.xls $dir2/$array[0]/$array[0]\_$ele");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/ccs.fasta $dir2/$array[0]/$array[0]\_$ele/read_of_insert.fasta");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/isoseq*fl*fasta $dir2/$array[0]/$array[0]\_$ele");
system("cp $dir1/$array[0]/$array[0]\_$ele/*/*consensus_isoforms* $dir2/$array[0]/$array[0]\_$ele");
}
}
###json to stat.xls 
$cDNA=join (",",@cDNAs);
system("python $Bin/IsoSeq/iso_stat_sum.py --indir $dir1/ --sample $array[0] --cDNA $cDNA --outdir $dir2/$array[0]");
system("$Bin/IsoSeq/merge_roi_mulsamples_data --indir $dir2/ -name roi_stat.xls --sample $array[0] --outfile $dir2/all_roi_stat.xls");
system("$Bin/IsoSeq/merge_roi_mulsamples_data --indir $dir2/ -name classify_stat.xls --sample $array[0] --outfile $dir2/all_classify_stat.xls");
system("$Bin/IsoSeq/merge_roi_mulsamples_data --indir $dir2/ -name cluster_stat.xls --sample $array[0] --outfile $dir2/all_cluster_stat.xls");
##

system("cat $dir2/$array[0]/*/roi_length.xls >$dir2/$array[0]/$array[0]\_roi_length.xls");
system("Rscript $Bin/IsoSeq/plot_length.R $dir2/$array[0]/$array[0]\_roi_length.xls $dir2/$array[0]/$array[0]\_roi_length_density");
system("cat $dir2/$array[0]/*/flnc_length.xls >$dir2/$array[0]/$array[0]\_flnc_length.xls");
system("Rscript $Bin/IsoSeq/plot_length.R $dir2/$array[0]/$array[0]\_flnc_length.xls $dir2/$array[0]/$array[0]\_flnc_length_density");
system("cat $dir2/$array[0]/*/consensus_length.xls >$dir2/$array[0]/$array[0]\_consensus_length.xls");
system("Rscript $Bin/IsoSeq/plot_length.R $dir2/$array[0]/$array[0]\_consensus_length.xls $dir2/$array[0]/$array[0]\_consensus_length_density");
# system("rm -r $dir2/$array[0]/*/*.json");
system("rm -r $dir2/$array[0]/*/*length.xls");
system("rm -r $dir2/$array[0]/*length.xls");
#####################################################################################
#                              collect the output of LSC                            #
#####################################################################################
if (-e "$path/LSC/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/LSC/$array[0]");
system("cp  $path/LSC/$array[0]/cat_lq_consensus_isoforms.fasta $path/PacBio_RNA_denovo_report/result/LSC/$array[0]/$array[0]\_all_low_qv_consensus_isoforms.fasta");
system("cp  $path/LSC/$array[0]/full_LR.fa $path/PacBio_RNA_denovo_report/result/LSC/$array[0]/$array[0]_LSC.fasta");
}
#####################################################################################
#                       collect the output of CD-hit-est                            #
#####################################################################################
if (-e "$path/CD-hit-est/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/CD-hit-est/$array[0]");
system("cp  $path/CD-hit-est/$array[0]/cat_hq_consensus_isoforms.fasta $path/PacBio_RNA_denovo_report/result/CD-hit-est/$array[0]/$array[0]\_all_high_qv_consensus_isoforms.fasta");
system("cp  $path/CD-hit-est/$array[0]/UniIso.fa $path/PacBio_RNA_denovo_report/result/CD-hit-est/$array[0]/CD-hit-est.fasta");
}
#####################################################################################
#                       collect the output of Transdecoder                          #
#####################################################################################
if (-e "$path/Transdecoder/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/Transdecoder/$array[0]");
system("cp  $path/Transdecoder/$array[0]/UniIso.fa.transdecoder.gff3 $path/PacBio_RNA_denovo_report/result/Transdecoder/$array[0]/transdecoder.gff3");
system("cp  $path/Transdecoder/$array[0]/UniIso.fa.transdecoder.cds.txt $path/PacBio_RNA_denovo_report/result/Transdecoder/$array[0]/transdecoder.cds.fa");
system("cp  $path/Transdecoder/$array[0]/UniIso.fa.transdecoder.pep.txt $path/PacBio_RNA_denovo_report/result/Transdecoder/$array[0]/transdecoder.pep.fa");
}
#####################################################################################
#                       collect the output of OrthoMCL                              #
#####################################################################################
if (-e "$path/OrthoMCL/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/OrthoMCL/$array[0]");
system("cp  $path/OrthoMCL/$array[0]/$prefix.group.xls $path/PacBio_RNA_denovo_report/result/OrthoMCL/$array[0]/group.xls");
}
#####################################################################################
#                       collect the output of Annotation                            #
#####################################################################################
if (-e "$path/Annotation/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]");
system("mkdir $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/GO");
system("mkdir $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KEGG");
system("mkdir $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KOG");
system("mkdir $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/NT_NR_SW");
system("cp $path/Annotation/$array[0]/Result/Function_Annotation.stat.xls $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]");
system("cp $path/Annotation/$array[0]/Result/Integrated_Function.annotation.xls $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]");
system("cp $path/Annotation/$array[0]/Result/UniIso.fa.Anno_Venn.png $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/Anno_Venn.png");
system("cp $path/Annotation/$array[0]/Result/GO_classification* $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/GO");
system("cp $path/Annotation/$array[0]/Result/UniIso.fa.GO.anno.xls $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/GO/GO.anno.xls");
system("cp $path/Annotation/$array[0]/Result/KEGG_classification_count.txt $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KEGG");
system("sed -i '1d' $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KEGG/KEGG_classification_count.txt");
system("cp $path/Annotation/$array[0]/Result/KEGG_classification.pdf $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KEGG");
system("cp $path/Annotation/$array[0]/Result/KEGG_classification.png $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KEGG");
system("cp $path/Annotation/$array[0]/Result/UniIso.fa.KEGG.gene2ko.xls $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KEGG/KEGG.gene2ko.xls");
system("cp $path/Annotation/$array[0]/Result/*KOG*class* $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/KOG");
system("cp $path/Annotation/$array[0]/Result/UniIso.fa.NR.anno.xls $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/NT_NR_SW/NR.anno.xls");
system("cp $path/Annotation/$array[0]/Result/UniIso.fa.NT.anno.xls $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/NT_NR_SW/NT.anno.xls");
system("cp $path/Annotation/$array[0]/Result/UniIso.fa.Swissprot.anno.xls $path/PacBio_RNA_denovo_report/result/Annotation/$array[0]/NT_NR_SW/Swissprot.anno.xls");
}
#####################################################################################
#                       collect the output of SSR                                   #
#####################################################################################
if (-e "$path/SSR/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/SSR/$array[0]");
system("cp  $path/SSR/$array[0]/SSR_Primer/UniIso.fa.misa $path/PacBio_RNA_denovo_report/result/SSR/$array[0]");
system("cp  $path/SSR/$array[0]/SSR_Primer/$prefix.results.xls $path/PacBio_RNA_denovo_report/result/SSR/$array[0]");
system("cp  $path/SSR/$array[0]/SSR_Primer/*.txt $path/PacBio_RNA_denovo_report/result/SSR/$array[0]");
system("cp  $path/SSR/$array[0]/SSR_Primer/ $path/PacBio_RNA_denovo_report/result/SSR/$array[0]");
system("cp  $path/SSR/$array[0]/SSR_Primer/Transcripts.fa.SSR_motifs_distribution.* $path/PacBio_RNA_denovo_report/result/SSR/$array[0]");
}
#####################################################################################
#                       collect the output of LncRNA                                #
#####################################################################################
if (-e "$path/LncRNA/$array[0]"){
system("mkdir $path/PacBio_RNA_denovo_report/result/LncRNA/$array[0]");
system("cp  $path/LncRNA/$array[0]/UniIso_ncrna.fa $path/PacBio_RNA_denovo_report/result/LncRNA/$array[0]/$array[0]_lncrna.fa");
system("cp  $path/LncRNA/$array[0]/lncRNA_length* $path/PacBio_RNA_denovo_report/result/LncRNA/$array[0]");
}
}
close IN;
#####################################################################################
#                       collect the output of RefIso                                #
#####################################################################################
if (-e "$path/RefIso"){
system("cp $path/RefIso/RefIso.fa $path/PacBio_RNA_denovo_report/result/RefIso");
if (-e "$path/RefIso/Multiple_isoforms.fasta"){
system("cp $path/RefIso/Multiple_isoforms.fasta $path/PacBio_RNA_denovo_report/result/RefIso/pre_RefIso.fa");}
}
#####################################################################################
#                       collect the output of expression evaluation                 #
#####################################################################################
if (-e "$path/RSEM"){
open SH,"$path/RSEM/sample_list.txt";
my $line=<SH>;
my $combine="";
chomp $line;
my @sample=split /,/,$line;
foreach my $ele(@sample){$combine=$combine." ".$path."/RSEM/".$ele."/".$ele.".Readcount_FPKM.xls";}
system("$Bin/expression/merge_FPKM_single_table $combine >$path/RSEM/merged.fpkm.xls");
system("$Bin/expression/merge_readcounts_single_table $combine >$path/RSEM/merged.readcount.xls");
system("$Bin/expression/calrowmeans -rpkm $path/RSEM/merged.fpkm.xls -group $line -groupname $line -out-rowmeans $path/RSEM/rowmeans_fpkm.xls");
system("$Bin/expression/plot.FPKM_v3.1 -r $path/RSEM/rowmeans_fpkm.xls -fpkm $path/RSEM/merged.fpkm.xls -rc $path/RSEM/merged.readcount.xls -len $path/RSEM/$sample[0]/geneINFO -output $path/RSEM");
my $dir = "$path/RSEM";
my @files=();
sub wanted {
    if ( -d $File::Find::dir ) {
        if ( $File::Find::dir =~ /^$path\/RSEM\/([^\/]+)$/ )#give the standard
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
system("mkdir  $path/PacBio_RNA_denovo_report/result/Expression_evaluation/$newpath[1]");
system("cp  $path/RSEM/merged.fpkm.xls $path/PacBio_RNA_denovo_report/result/Expression_evaluation/");
system("cp  $path/RSEM/FPKM_boxplot* $path/PacBio_RNA_denovo_report/result/Expression_evaluation/");
system("cp  $path/RSEM/FPKM_density_distribution* $path/PacBio_RNA_denovo_report/result/Expression_evaluation/");
system("cp  $path/RSEM/*/*TranscriptExp.xls $path/PacBio_RNA_denovo_report/result/Expression_evaluation/$newpath[1]");
system("cp  $path/RSEM/*/*.mapped.stat.xls $path/PacBio_RNA_denovo_report/result/Expression_evaluation/$newpath[1]");
}
system("cat $path/RSEM/*/*.mapped.stat.2.xls | sed '1i Sample	Total Reads	Mapped Reads(ratio)	Uniq mapped Reads(ratio)	Multi mapped Reads(ratio)' > $path/PacBio_RNA_denovo_report/result/Expression_evaluation/Total_mapped.stat.xls");
}
#####################################################################################
#                       collect the output of differential expression               #
#####################################################################################
if (-e "$path/edgeR"){
	system("cp $path/edgeR/edgeR.*.dir/*edgeR.DE_results_diff_exp.xls $path/PacBio_RNA_denovo_report/result/Differential_expression");
#	system("cp $path/edgeR/edgeR.*.dir/*edgeR.DE_results.MA* $path/PacBio_RNA_denovo_report/result/Differential_expression");
	system("cp $path/edgeR/edgeR.*.dir/*edgeR.DE_results.Volcano* $path/PacBio_RNA_denovo_report/result/Differential_expression");
	system("cp $path/edgeR/edgeR.*.dir/*genes_vs_samples_heatmap.png $path/PacBio_RNA_denovo_report/result/Differential_expression");
	system("cp $path/edgeR/edgeR.*.dir/*sample_cor_matrix.png $path/PacBio_RNA_denovo_report/result/Differential_expression");
	system("cp -r $path/edgeR/edgeR.*.dir/GO $path/PacBio_RNA_denovo_report/result/Differential_expression");
	system("cp -r $path/edgeR/edgeR.*.dir/KEGG/ $path/PacBio_RNA_denovo_report/result/Differential_expression");
	system("rm -r $path/PacBio_RNA_denovo_report/result/Differential_expression/KEGG/*/*.annotate");
	system("rm -r $path/PacBio_RNA_denovo_report/result/Differential_expression/KEGG/*/*.blast");
	system("rm -r $path/PacBio_RNA_denovo_report/result/Differential_expression/KEGG/*/*.identify");
	system("rm -r $path/PacBio_RNA_denovo_report/result/Differential_expression/KEGG/*/*.Kegg.ko2gene.xls");
	system("rm -r $path/PacBio_RNA_denovo_report/result/Differential_expression/KEGG/*/*.tmp.xls");
}
#####################################################################################
#                               collect the output of AS                            #    
#####################################################################################

if (-e "$path/AS"){
	system("cp -r $path/AS/result $path/PacBio_RNA_denovo_report/result/AS");
}
if (-e "$path/AS/GO"){
	system("cp -r $path/AS/GO $path/PacBio_RNA_denovo_report/result/AS");
}
if (-e "$path/AS/KEGG"){
	system("cp -r $path/AS/KEGG $path/PacBio_RNA_denovo_report/result/AS");
	system("rm -r $path/PacBio_RNA_denovo_report/result/AS/KEGG/*.annotate");
	system("rm -r $path/PacBio_RNA_denovo_report/result/AS/KEGG/*.blast");
	system("rm -r $path/PacBio_RNA_denovo_report/result/AS/KEGG/*.identify");
	system("rm -r $path/PacBio_RNA_denovo_report/result/AS/KEGG/*.Kegg.ko2gene.xls");
	system("rm -r $path/PacBio_RNA_denovo_report/result/AS/KEGG/*.tmp.xls");
}
#####################################################################################
#                                  produce the web report                           #
#####################################################################################
	open SH, ">$path/PacBio_RNA_denovo_report/web.sh";
	print SH "/share/public/software/node-v6.11.1/bin/node $Bin/jade2html/pac_rna_denovo.js ./";
	close SH;
	system("cd $path/PacBio_RNA_denovo_report/ && qsub -cwd web.sh");
####################################################################################
#                                 END                                              #
####################################################################################                                 

