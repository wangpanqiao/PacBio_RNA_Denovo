#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use FindBin qw($Bin $Script);
use File::Basename;
use lib ("/share/public/software/trinity/trinityrnaseq_r20140717/PerlLib/");
use Fasta_reader;
use Data::Dumper;


my $usage = <<__EOUSAGE__;


#################################################################################################
#
#  Required:
#
#  --matrix <string>               matrix of raw read counts (not normalized!)
#
#  --method <string>               edgeR|DESeq|DESeq2   (DESeq(2?) only supported here w/ bio replicates)
#                                   
#
#  Optional:
#
#  --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1
#                                        cond_A    cond_A_rep2
#                                        cond_B    cond_B_rep1
#                                        cond_B    cond_B_rep2
#  --FDR     <float>               default: 0.05  fdr threshold for differential expressed gene corrected pvalue
#  --log2FC  <int>                 default: 1     log2FC threshold (log2 of gene express foldchange)
#
#  General options:
#
#  --min_rowSum_counts <int>       default: 2  (only those rows of matrix meeting requirement will be tested)
#
#  --output|o                      name of directory to place outputs (default: \$method.\$pid.dir)
#
#  --reference_sample <string>     name of a sample to which all other samples should be compared.
#                                   (default is doing all pairwise-comparisons among samples)
#
#  --contrasts <string>            file (tab-delimited) containing the pairs of sample comparisons to perform.
#                                  ex. 
#                                       cond_A    cond_B
#                                       cond_Y    cond_Z
#
#
###############################################################################################
#
#  ## EdgeR-related parameters
#  ## (no biological replicates)
#
#  --dispersion <float>            edgeR dispersion value (default: 0.1)   set to 0 for poisson (sometimes breaks...)
#  --edgeR_fdr_method  <string>    edgeR false discovery control method (default: BH) ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#                                  DESeq only support "BH" method for FDR control
#
#  http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
#
###############################################################################################
#
#  ## DE-Seq related parameters
#
#  --DESEQ_method <string>         "pooled", "pooled-CR", "per-condition", "blind"
#  --DESEQ_sharingMode <string>    "maximum", "fit-only", "gene-est-only"
#  --DESEQ_fitType <string>        fitType = c("parametric", "local")
#
#  ## (no biological replicates)
#        note: FIXED as: method=blind, sharingMode=fit-only
#  ## (with biological replicates)
#        note: Default: method=pooled-CR, sharingMode=fit-only, fitType=parametric
#
#  http://www.bioconductor.org/packages/release/bioc/html/DESeq.html
#
################################################################################################



__EOUSAGE__


    ;


my ($matrix_file,$method,$samples_file,$min_rowSum_counts,$help_flag,$output_dir,$dispersion);
my ($edgeR_fdr_method,$fdr,$log2FC,$contrasts_file,$reference_sample,$DESEQ_method,$DESEQ_sharingMode,$DESEQ_fitType);
&GetOptions	( 'h' => \$help_flag,
			  'matrix=s' => \$matrix_file,
			  'method=s' => \$method,
			  'samples_file=s' => \$samples_file,
			  'FDR=f' => \$fdr,
			  'log2FC=i' => \$log2FC,
			  'output|o=s' => \$output_dir,
			  'min_rowSum_counts=i' => \$min_rowSum_counts,

			  'dispersion=f' => \$dispersion,
			  'edgeR_fdr_method=s' => \$edgeR_fdr_method,
	
			  'reference_sample=s' => \$reference_sample,
			  'contrasts=s' => \$contrasts_file,
			  

			  'DESEQ_method=s' => \$DESEQ_method,
			  'DESEQ_sharingMode=s' => \$DESEQ_sharingMode,
			  'DESEQ_fitType=s' => \$DESEQ_fitType,


	);
	
$min_rowSum_counts ||= 2;
$dispersion ||= 0.1;
$edgeR_fdr_method ||= "BH";
$fdr ||= "0.05";
$log2FC ||= "1";
$DESEQ_method ||= "pooled-CR";
$DESEQ_sharingMode ||= "fit-only";
$DESEQ_fitType ||= "parametric";

if ($help_flag) {
	die $usage;
}


unless ($matrix_file && $method) { 
	die $usage;
}

if ($matrix_file =~ /fpkm/i) {
	die "Error, be sure you're using a matrix file that corresponds to raw counts, and not FPKM values.\n" . "If this is correct, then please rename your file, and remove fpkm from the name.\n\n";
}


unless ($method =~ /^(edgeR|DESeq2?|GLM)$/) {
	die "Error, do not recognize method: [$method], only edgeR or DESeq currently.";
}

main: {


	my $workdir = cwd();
	
	
	my %sample_name_to_column = &get_sample_name_to_column_index($matrix_file);
	
	my %samples;
	if ($samples_file) {
		unless ($samples_file =~ /^\//) {
			$samples_file = cwd() . "/$samples_file";
		}
		
		%samples = &parse_sample_info($samples_file);
	}
	else {
		# no replicates, so	assign each	sample to itself as	a single replicate
		foreach my $sample_name (keys %sample_name_to_column) {
			$samples{$sample_name} = [$sample_name];
		}
	}

	print Dumper(\%samples);
		
	if ($matrix_file !~ /^\//) {
		## make full path
		$matrix_file = cwd() . "/$matrix_file";
	}
		
	unless ($output_dir) {
		$output_dir = "$method.$$.dir";
	}

	`rm -rf $output_dir` if (-d $output_dir);
	mkdir($output_dir) or die "Error, cannot mkdir $output_dir";
	chdir $output_dir or die "Error, cannot cd to $output_dir";
	

	my @sample_names = keys %samples;


	if ($method eq "GLM") {
		unless ($samples_file) { 
			die "Error, need samples file for GLM";
		}
		## samples file here requires a different format:
		# replicate (tab) attrA [(tab) attrB, ...]

		&run_GLM($matrix_file, \%samples, \%sample_name_to_column);
	}
	else {
		# edgeR or DESeq pairwise comparison between samples:
		
		my @DE_contrasts;
		
		if ($reference_sample) {
			
			my @other_samples = grep { $_ ne $reference_sample} @sample_names;

			unless (@other_samples) {
				die "Error, couldn't extract non-reference samples from list: @sample_names";
			}

			foreach my $other_sample (@other_samples) {
				push (@DE_contrasts, [$reference_sample, $other_sample]);
			}
			
		}
		elsif ($contrasts_file) {
			
			unless ($contrasts_file =~ /^\//) {
				$contrasts_file = "$workdir/$contrasts_file";
			}
			
			open (my $fh, $contrasts_file) or die "Error, cannot open file $contrasts_file";
			while (<$fh>) {
				chomp;
				unless (/\w/) { next; }
				if (/^\#/) { next; }
				my ($sampleA, $sampleB) = split(/\t/);
				unless ($sampleA && $sampleB) {
					die "Error, didn't read a pair of tab-delimited samples from $contrasts_file, line: $_";
				}
				push (@DE_contrasts, [$sampleA, $sampleB]);
			}
			close $fh;
		}
		else {
			## performing all pairwise comparisons:
			
			@sample_names = sort @sample_names;
			for (my $i = 0; $i < $#sample_names; $i++) {
				for (my $j = $i + 1; $j <= $#sample_names; $j++) {

					push (@DE_contrasts, [$sample_names[$i], $sample_names[$j]]);
				}
			}
		}
		
		print STDERR "Contrasts to perform are: " . Dumper(\@DE_contrasts);
		
		foreach my $DE_contrast (@DE_contrasts) {

			my ($sample_a, $sample_b) = @$DE_contrast;
			
			if ($method	eq "edgeR") {
				&run_edgeR_sample_pair($matrix_file, \%samples, \%sample_name_to_column, $sample_a, $sample_b);
				
			}
			elsif ($method eq "DESeq") {
				&run_DESeq_sample_pair($matrix_file, \%samples,	\%sample_name_to_column, $sample_a, $sample_b);
			}
			elsif ($method eq "DESeq2") {
				&run_DESeq2_sample_pair($matrix_file, \%samples, \%sample_name_to_column, $sample_a, $sample_b);
			}

		}
	}
		

	exit(0);
}

####
sub parse_sample_info {
	my ($sample_file) = @_;

	my %samples;

	open (my $fh, $sample_file) or die $!;
	while (<$fh>) {
		unless (/\w/) { next; }
		if (/^\#/) { next; } # allow comments
		chomp;
		s/^\s+//; # trim any leading ws
		my @x = split(/\s+/); #	now	ws instead of just tabs
		if (scalar @x < 2) { next; }
		my ($sample_name, $replicate_name, @rest) = @x;
		
		#$sample_name =~ s/^\s|\s+$//g;
		#$replicate_name =~ s/^\s|\s+$//g;
		
		push (@{$samples{$sample_name}}, $replicate_name);
	}
	close $fh;

	return(%samples);
}

####
sub	get_sample_name_to_column_index	{
	my ($matrix_file) = @_;

	my %column_index;

	open (my $fh, $matrix_file) or die "Error, cannot open file $matrix_file";
	my $header_line = <$fh>;

	$header_line =~ s/^\#//; # remove comment field.
	$header_line =~ s/^\s+|\s+$//g;
	my @samples = split(/\t/, $header_line);

	{ # check for disconnect between header	line and data lines
		my $next_line = <$fh>;
		my @x = split(/\t/, $next_line);
		print STDERR "Got " . scalar(@samples) . " samples, and got: " . scalar(@x) . " data fields.\n";
		print STDERR "Header: $header_line\nNext: $next_line\n";
		
		if (scalar(@x) == scalar(@samples)) {
			# problem... shift headers over, no	need for gene column heading
			shift @samples;
			print STDERR "-shifting	sample indices over.\n";
		}
	}
	close $fh;
			
	
	my $counter = 0;
	foreach my $sample (@samples) {
		$counter++;
		
		$sample	=~ s/\.(isoforms|genes)\.results$//; 
		
		$column_index{$sample} = $counter;
	}

	use Data::Dumper;
	print STDERR Dumper(\%column_index);
	

	return(%column_index);
	
}


####
sub	run_edgeR_sample_pair {
	my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B) = @_;
		 
	my $output_prefix = join("_vs_", ($sample_A, $sample_B));
		
	my $Rscript_name = "$output_prefix.$sample_A.vs.$sample_B.EdgeR.Rscript";

	my @reps_A = @{$samples_href->{$sample_A}};
	my @reps_B = @{$samples_href->{$sample_B}};

	my $num_rep_A = scalar(@reps_A);
	my $num_rep_B = scalar(@reps_B);
	
	my @rep_column_indices;
	foreach my $rep_name (@reps_A, @reps_B) {
		my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot	determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
		push (@rep_column_indices, $column_index);
	}
		

	## write R-script to run edgeR
	open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
	
#	print $ofh ".libPaths(\"/share/work1/staff/mengfei/R_install_packges/\")\n";
	print $ofh "library(edgeR)\n";

	print $ofh "\n";
	
	print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
	print $ofh "col_ordering = c(" . join(",", @rep_column_indices) . ")\n";
	print $ofh "rnaseqMatrix = data[,col_ordering]\n";
	print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
	print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=$min_rowSum_counts,]\n";
	print $ofh "conditions = factor(c(rep(\"$sample_A\", $num_rep_A), rep(\"$sample_B\", $num_rep_B)))\n";
	print $ofh "\n";
	print $ofh "exp_study = DGEList(counts=rnaseqMatrix, group=conditions)\n";
	print $ofh "exp_study = calcNormFactors(exp_study)\n";
	
	if ($num_rep_A > 1 && $num_rep_B > 1) {
		print $ofh "exp_study = estimateCommonDisp(exp_study)\n";
		print $ofh "exp_study = estimateTagwiseDisp(exp_study)\n";
		print $ofh "et = exactTest(exp_study)\n";
	}
	else {
		print $ofh "et = exactTest(exp_study, dispersion=$dispersion)\n";
	}
	print $ofh "tTags = topTags(et,n=NULL,adjust.method=\"$edgeR_fdr_method\",sort.by=\"none\")\n";
	print $ofh "plot_data = tTags\n";
#	if ($num_rep_A > 1 && $num_rep_B > 1) {
#		print $ofh "pseudo_counts = exp_study\$pseudo.counts\n";
#		print $ofh "tTags = cbind(pseudo_counts, as.data.frame(tTags))\n";
#	}
#	else {
#		print $ofh "raw_counts = exp_study\$counts\n";
#		print $ofh "tTags = cbind(raw_counts, as.data.frame(tTags))\n";
#	}

	print $ofh "write.table(tTags, file=\'$output_prefix.edgeR.DE_results\', sep='\t', quote=F, row.names=T)\n";

	## generate	MA and Volcano plots
	print $ofh	"source(\"$Bin/R/rnaseq_plot_funcs.R\")\n";
	print $ofh	"result_table = plot_data\$table\n";
	print $ofh	"png(\"$output_prefix.edgeR.DE_results.MA.png\")\n";
	print $ofh	"plot_MA(result_table\$logCPM, result_table\$logFC, result_table\$FDR, $fdr, $log2FC, pch='.')\n";
	print $ofh	"png(\"$output_prefix.edgeR.DE_results.Volcano.png\")\n";
	print $ofh	"plot_Volcano(result_table\$logFC, result_table\$FDR, $fdr, $log2FC, pch='.')\n";
	print $ofh	"pdf(\"$output_prefix.edgeR.DE_results.MA.pdf\")\n";
	print $ofh	"plot_MA(result_table\$logCPM, result_table\$logFC, result_table\$FDR, $fdr, $log2FC, pch='.')\n";
	print $ofh	"pdf(\"$output_prefix.edgeR.DE_results.Volcano.pdf\")\n";
	print $ofh	"plot_Volcano(result_table\$logFC, result_table\$FDR, $fdr, $log2FC, pch='.')\n";
	print $ofh	"dev.off()\n";
	
	close $ofh;

	## Run R-script
	my $cmd = "R --vanilla -q < $Rscript_name";


	eval {
		&process_cmd($cmd);
	};
	if ($@) {
		print STDERR "$@\n\n";
		print STDERR "\n\nWARNING: This EdgeR comparison failed...\n\n";
		## if this is due to data paucity, such	as in small	sample data	sets, then ignore for now.
	}
	

	return;
}

sub	run_DESeq_sample_pair {
	my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B) = @_;
	
	my $output_prefix = join("_vs_", ($sample_A, $sample_B));
	
	my $Rscript_name = "$output_prefix.DESeq.Rscript";
	
	my @reps_A = @{$samples_href->{$sample_A}};
	my @reps_B = @{$samples_href->{$sample_B}};

	my $num_rep_A = scalar(@reps_A);
	my $num_rep_B = scalar(@reps_B);
	
	
	my @rep_column_indices;
	foreach	my $rep_name (@reps_A, @reps_B) {
		my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
		push (@rep_column_indices, $column_index);
	}
	

	## write R-script to run edgeR
	open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
	
	print $ofh ".libPaths(\"/share/work1/staff/mengfei/R_install_packges/\")\n";
	print $ofh "library(DESeq)\n";
	print $ofh "\n";

	print $ofh "data = read.table(\"$matrix_file\",	header=T, row.names=1, com='')\n";
	print $ofh "col_ordering = c(" . join(",", @rep_column_indices)	. ")\n";
	print $ofh "rnaseqMatrix = data[,col_ordering]\n";
	print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
	print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=$min_rowSum_counts,]\n";
	print $ofh "conditions = factor(c(rep(\"$sample_A\", $num_rep_A), rep(\"$sample_B\", $num_rep_B)))\n";
	print $ofh "\n";
	print $ofh "exp_study = newCountDataSet(rnaseqMatrix, conditions)\n";
	print $ofh "exp_study = estimateSizeFactors(exp_study)\n";
	#print $ofh	"sizeFactors(exp_study)\n";
	#print $ofh	"exp_study = estimateVarianceFunctions(exp_study)\n";
	
	if ($num_rep_A == 1 && $num_rep_B == 1) {
		
		print STDERR "\n\n** Note, no replicates, setting method='blind', sharingMode='fit-only'\n\n";
		
		$DESEQ_method = "blind";
		$DESEQ_sharingMode = "fit-only";
		
	}

	# got bio replicates
	my $est_disp_cmd = "exp_study = estimateDispersions(exp_study";
	
	if ($DESEQ_method) {
		$est_disp_cmd .= ", method=\"$DESEQ_method\"";
	}
	
	if ($DESEQ_sharingMode) {
		$est_disp_cmd .= ", sharingMode=\"$DESEQ_sharingMode\"";
	}
	
	if ($DESEQ_fitType) {
		$est_disp_cmd .= ", fitType=\"$DESEQ_fitType\"";
	}
	
	$est_disp_cmd .= ")\n";
	
	print $ofh $est_disp_cmd;
	
	
	#print $ofh	"str(fitInfo(exp_study))\n";
	#print $ofh	"plotDispEsts(exp_study)\n";
	print $ofh "\n";
	print $ofh "res = nbinomTest(exp_study, \"$sample_A\", \"$sample_B\",)\n";
	print $ofh "\n";
## output results
	print $ofh "write.table(res[order(res\$pval),], file=\'$output_prefix.DESeq.DE_results\', sep='\t', quote=FALSE, row.names=FALSE)\n";
	
	## generate MA and Volcano plots
	print $ofh "source(\"$Bin/R/rnaseq_plot_funcs.R\")\n";
#	print $ofh "pdf(\"$output_prefix.DESeq.DE_results.MA_n_Volcano.pdf\")\n";
#	print $ofh "plot_MA_and_Volcano(log2(res\$baseMean+1),	res\$log2FoldChange, res\$padj)\n";

	print $ofh "png(\"$output_prefix.DESeq.DE_results.MA.png\")\n";
	print $ofh "plot_MA(log2(res\$baseMean+1), res\$log2FoldChange, res\$padj, $fdr, $log2FC, pch='.')\n";
	print $ofh "png(\"$output_prefix.DESeq.DE_results.Volcano.png\")\n";
	print $ofh "plot_Volcano(res\$log2FoldChange, res\$padj, $fdr, $log2FC, pch='.')\n";

	print $ofh "pdf(\"$output_prefix.DESeq.DE_results.MA.pdf\")\n";
	print $ofh "plot_MA(log2(res\$baseMean+1), res\$log2FoldChange, res\$padj, $fdr, $log2FC, pch='.')\n";
	print $ofh "pdf(\"$output_prefix.DESeq.DE_results.Volcano.pdf\")\n";
	print $ofh "plot_Volcano(res\$log2FoldChange, res\$padj, $fdr, $log2FC, pch='.')\n";
	print $ofh "dev.off()\n";

	close $ofh;
	


	## Run R-script
	my $cmd = "R --vanilla -q < $Rscript_name";
	&process_cmd($cmd);
	&DESeq_result_format_ ("$output_prefix.DESeq.DE_results",$fdr,$log2FC,$output_prefix,$sample_A,$sample_B);

	return;
}

sub run_DESeq2_sample_pair {
	my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B) = @_;
		 
	my $output_prefix =	join("_vs_", ($sample_A, $sample_B));
		
	my $Rscript_name = "$output_prefix.DESeq2.Rscript";
	
	my @reps_A = @{$samples_href->{$sample_A}};
	my @reps_B = @{$samples_href->{$sample_B}};

	my $num_rep_A = scalar(@reps_A);
	my $num_rep_B = scalar(@reps_B);
	

	if ($num_rep_A < 2 || $num_rep_B < 2) {
		print STDERR "DESeq2 only supported here with biological replicates for each condition. Skipping: $sample_A vs. $sample_B *** \n\n";
		return;
	}
	
	my @rep_column_indices;
	foreach	my $rep_name (@reps_A, @reps_B) {
		my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot	determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
		push (@rep_column_indices, $column_index);
	}
	

	## write R-script to run edgeR
	open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
	
##	print $ofh ".libPaths(\"/share/work1/staff/mengfei/R_install_packges/\")\n";
	print $ofh "library(DESeq2)\n";
	print $ofh "\n";

	print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
	print $ofh "col_ordering = c(" . join(",", @rep_column_indices) . ")\n";
	print $ofh "rnaseqMatrix = data[,col_ordering]\n";
	print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
	print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=$min_rowSum_counts,]\n";
	print $ofh "conditions = data.frame(conditions=factor(c(rep(\"$sample_A\", $num_rep_A),	rep(\"$sample_B\", $num_rep_B))))\n";
	print $ofh "rownames(conditions) = colnames(rnaseqMatrix)\n";
	print $ofh "ddsFullCountTable <- DESeqDataSetFromMatrix(\n"
			 . "    countData = rnaseqMatrix,\n"
			 . "    colData = conditions,\n"
			 . "    design = ~ conditions)\n";
	print $ofh "dds = DESeq(ddsFullCountTable)\n";
	print $ofh "res = results(dds)\n";
	print $ofh "plot_data = res\n";
	print $ofh "baseMeanA = rowMeans(counts(dds,normalized=TRUE)[,colData(dds)\$condition == \"$sample_A\"])\n";
	print $ofh "baseMeanB = rowMeans(counts(dds,normalized=TRUE)[,colData(dds)\$condition == \"$sample_B\"])\n";
	print $ofh "res = cbind(baseMeanA, baseMeanB, as.data.frame(res))\n";

	print $ofh "\n";
	## output results
	print $ofh "write.table(as.data.frame(res[order(res\$pvalue),]), file=\'$output_prefix.DESeq2.DE_results\', sep='\t', quote=FALSE, row.names=T)\n";
	
#	##	generate MA and Volcano plots
	print $ofh "source(\"$Bin/R/rnaseq_plot_funcs.R\")\n";
#	print $ofh "pdf(\"$output_prefix.DESeq2.DE_results.MA_n_Volcano.pdf\")\n";
#	print $ofh "plot_MA_and_Volcano(log2(res\$baseMean+1), res\$log2FoldChange, res\$padj)\n";

	print $ofh "png(\"$output_prefix.DESeq2.DE_results.MA.png\")\n";
	print $ofh "plot_MA(log2(plot_data\$baseMean+1), plot_data\$log2FoldChange, plot_data\$padj, $fdr, $log2FC, pch='.')\n";
	print $ofh "png(\"$output_prefix.DESeq2.DE_results.Volcano.png\")\n";
	print $ofh "plot_Volcano(plot_data\$log2FoldChange, plot_data\$padj, $fdr, $log2FC, pch='.')\n";

	print $ofh "pdf(\"$output_prefix.DESeq2.DE_results.MA.pdf\")\n";
	print $ofh "plot_MA(log2(plot_data\$baseMean+1), plot_data\$log2FoldChange, plot_data\$padj, $fdr, $log2FC, pch='.')\n";
	print $ofh "pdf(\"$output_prefix.DESeq2.DE_results.Volcano.pdf\")\n";
	print $ofh "plot_Volcano(plot_data\$log2FoldChange, plot_data\$padj, $fdr, $log2FC, pch='.')\n";
	print $ofh "dev.off()\n";
	
	
	close $ofh;
	
	## Run R-script
	my $cmd = "R --vanilla -q < $Rscript_name";
	&process_cmd($cmd);
	&DESeq2_result_format_ ("$output_prefix.DESeq2.DE_results",$fdr,$log2FC,$output_prefix,$sample_A,$sample_B);

	return;
}




####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret) ";
	}

	return;
}

####
sub	run_GLM {
	my ($matrix_file, $samples_href, $sample_name_to_column_index_href) = @_;
	

	my $output_prefix = basename($matrix_file);

			 
	my $Rscript_name = "$output_prefix.GLM.Rscript";
	
	## write R-script to run edgeR
	open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
	
	print $ofh ".libPaths(\"/share/work1/staff/mengfei/R_install_packges/\")";
	print $ofh "library(edgeR)\n";

	print $ofh "\n";
	
	print $ofh "design_matrix = read.table(\"$samples_file\", header=T, row.names=1)\n";
	print $ofh "groups = factor(apply(design_matrix, 1, paste, collapse='.'))\n";
	print $ofh "design_matrix = cbind(design_matrix, groups=groups)\n";
	
	print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
	print $ofh "rnaseqMatrix = round(data)\n";
	print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=$min_rowSum_counts,]\n";

	print $ofh "design = model.matrix(~0+groups)\n";
	print $ofh "colnames(design) = levels(groups)\n";
	print $ofh "rownames(design) = rownames(design_matrix)\n";
	print $ofh "rnaseqMatrix = rnaseqMatrix[,rownames(design)] # ensure properly ordered according to design\n";
	print $ofh "exp_study = DGEList(counts=rnaseqMatrix, group=groups)\n";
	print $ofh "exp_study = estimateGLMCommonDisp(exp_study,design)\n";
	print $ofh "exp_study = estimateGLMTrendedDisp(exp_study, design)\n";
	print $ofh "exp_study = estimateGLMTagwiseDisp(exp_study, design)\n";
	print $ofh "fit = glmFit(exp_study, design)\n";
	print $ofh "## define your contrasts:\n";
	print $ofh "levels(groups) # examine the factor combinations\n";
	print $ofh "# contrast = makeContrasts((wt.T15-wt.T0)-(zcf15.T15-zcf15.T0), levels=design)\n";
	print $ofh "# lrt = glmLRT(fit, contrast=contrast)\n";
	print $ofh "# topTags(lrt, n=100)\n";
	
	
	close $ofh;

	
	return;
}


sub edgeR_result_format_ {
	my ($result,$fdr,$log2FC,$output) = @_;

	open (my $in, "$result") or die $!;
	my $header_line = <$in>;
	open (my $out1, ">$output.Significant_DE.xls") or die $!;
	open (my $out2, ">$output.All_DE.xls") or die $!;
	print $out1 "#geneID\t$header_line";
	print $out2 "#geneID\t$header_line";
	while (<$in>) {
		chomp;
		s/\s+$//;
		next if (/^$/);
		my @info = split/\t+/,$_;
		if ($info[-1] <= $fdr && abs($info[-4]) >= $log2FC) {
			print $out1 "$_\n";
		}
		print $out2 "$_\n";
	}
	close $in;
	close $out1;
	close $out2;
}

sub DESeq_result_format_ {
	my ($result,$fdr,$log2FC,$output,$samA,$samB) = @_;

	open (my $in, "$result") or die $!;
	my $header_line = <$in>;
	open (my $out1, ">$output.Significant_DE.xls") or die $!;
	open (my $out2, ">$output.All_DE.xls") or die $!;
	print $out1 "#geneID\t$samA"."_norm_count\t$samB"."_norm_count\tlog2FC\tPValue\tFDR\n";
	print $out2 "#geneID\t$samA"."_norm_count\t$samB"."_norm_count\tlog2FC\tPValue\tFDR\n";
	while (<$in>) {
		chomp;
		s/\s+$//;
		next if (/^$/);
		my @info = split/\t+/,$_;
		if ($info[-1] <= $fdr && (abs($info[-3]) >= $log2FC || $info[-3] =~ /Inf/)) {
			print $out1 "$info[0]\t$info[2]\t$info[3]\t$info[5]\t$info[6]\t$info[7]\n";
		}
		print $out2 "$info[0]\t$info[2]\t$info[3]\t$info[5]\t$info[6]\t$info[7]\n";
	}
	close $in;
	close $out1;
	close $out2;
}

sub DESeq2_result_format_ {
	my ($result,$fdr,$log2FC,$output,$samA,$samB) = @_;

	open (my $in, "$result") or die $!;
	my $header_line = <$in>;
	open (my $out1, ">$output.Significant_DE.xls") or die $!;
	open (my $out2, ">$output.All_DE.xls") or die $!;
	print $out1 "#geneID\t$samA"."_norm_count\t$samB"."_norm_count\tlog2FC\tPValue\tFDR\n";
	print $out2 "#geneID\t$samA"."_norm_count\t$samB"."_norm_count\tlog2FC\tPValue\tFDR\n";
	while (<$in>) {
		chomp;
		s/\s+$//;
		next if (/^$/);
		my @info = split/\t+/,$_;
		print $out2 "$info[0]\t$info[1]\t$info[2]\t$info[4]\t$info[6]\t$info[7]\n";
		next if ($info[-1] =~ /NA/i);
		if ($info[-1] <= $fdr && abs($info[-4]) >= $log2FC) {
			print $out1 "$info[0]\t$info[1]\t$info[2]\t$info[4]\t$info[6]\t$info[7]\n";
		}
	}
	close $in;
	close $out1;
	close $out2;
}

