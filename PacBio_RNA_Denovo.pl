#!/usr/bin/perl

#  Copyright (c)   BerryGenomics-2016
#  Writer:         Xing Shilai <xingshilai3@163.com> <xingshilai547@berrygenomics.com>
#  Program Date:   2016.8.15
#  Modifier:       Xing Shilai <xingshilai3@163.com> <xingshilai547@berrygenomics.com>
#  Last Modified:  2016.9.5

use warnings;
use strict;
use Getopt::Long;
use Cwd qw (abs_path);
use FindBin qw ($Bin);
my $usage =
"
Usage:

  Options:
-project_name               <str>    The project name
-species                    <str>    The species name
-category                   <str>    the kegg category, animal|plant|fungi
-PacBio_samples_file        <file>   PacBio sample list with eight coloums and without header line, for example:\n sample name\tcDNA size\ts1_p0.1.bax.h5\ts1_p0.2.bax.h5\ts1_p0.3.bax.h5\tIllumina R1.clean.fastq.gz|NA\tIllumina R2.clean.fastq.gz|NA\tsample alias
-Iso_method		    <str>    Web|Command,default is Command
-corrected_reads_file       <file>   Illumina pair-end reads, for example:Pacbio sample\tA.R1.clean.fastq.gz\tA.R2.clean.fastq.gz\tB.R1.clean.fastq.gz\tetc
-RNAseq_samples_file        <file>   Illumina sample list file with three coloums and without header line, for example:\n sample name\tR1.clean.fastq.gz\tR2.clean.fastq.gz|NA
-group_file	            <file>   The case or control information for samples in RNAseq_samples_file 
-group_compare_file	    <file>   The comparision between two groups 
-parameters_file            <file>   Parameters_file, 
                                     for IsoSeq, 
                                     for LSC,
                                     for CD-hit-est, 
                                     for Annotation,
                                     for ORF transdecoder, 
                                     for OrthoMCL, 
                                     for SSR, 
                                     for LncRNA, 
                                     for Expression evaluation, 
                                     for Differential expression                                
-database_file              <file>   Database_file list for blast annotation
-DE_method	            <str>    edgeR|DESeq|DESeq2,default is edgeR
-prefix                     <str>    The prefix of out, (Species Name. etc)
-output_dir                 <str>    Output directionary, default is ./
-s|step                     <num>    Which step to run from,default is 1 [option: 1,2,3,4,5,6,7,8,9,10,11,12]
                                     1:IsoSeq,
                                     2:LSC,
				     3:CD-hit-est,
                                     4:Annotation,
                                     5:ORF transdecoder,
                                     6:OrthoMCL,
                                     7:SSR,
                                     8:LncRNA,
                                     9:Get reference isoform,
                                    10:Expression evaluation,
                                    11:Differential expression,
                                    12:Package result and web report
-cut                        <num>   Which step to end-running,default is 1 [option: 1,2,3,4,5,6,7,8,9,10,11,12]
-skip                       <num>   Which steps between step and cut to be ignored,sperated by comma
-h                                  Help
For example:
perl $0  --project_name  BFXXX  --species Mouse --category animal --PacBio_samples_file PacBio_samples_file.txt --RNAseq_samples_file RNAseq_samples_file.txt --group_file group_file.txt --group_compare_file group_compare_file.txt --parameters_file parameters_file.txt --database_file blast_database.txt --prefix Mouse --output_dir ./ --step 1 --cut 1    
";

my ($project_name, $species, $category, $PacBio_samples_file,$Iso_method, $corrected_reads_file, $RNAseq_samples_file, $group_file, $group_compare_file, $DE_method, $parameters_file, $database_file, $prefix, $output_dir, $step, $cut, $skip, $help);   

GetOptions(
                                    "project_name=s" => \$project_name,
	                            "species=s" => \$species,
	                            "category=s" => \$category,
	                            "PacBio_samples_file=s" => \$PacBio_samples_file,
	                            "Iso_method=s" => \$Iso_method,
	                            "corrected_reads_file=s" => \$corrected_reads_file,
	                            "RNAseq_samples_file=s" => \$RNAseq_samples_file,
	                            "group_file=s" => \$group_file,
	                            "group_compare_file=s" => \$group_compare_file,
	                            "DE_method=s" => \$DE_method,
	                            "parameters_file=s" => \$parameters_file,
	                            "database_file=s" => \$database_file,
	                            "prefix=s" => \$prefix,
	                            "output_dir=s" => \$output_dir,
	                            "s|step=i"=>\$step,
                                    "cut=i"=>\$cut,
                                    "skip=s"=>\$skip,
	                            "h|help"=>\$help
);
$output_dir ||= './';
$step ||= 1;
$cut ||= 1;
$skip||= 1;
$Iso_method ||= 'Command';
$DE_method ||= 'edgeR';
$output_dir = abs_path($output_dir);
$PacBio_samples_file = abs_path($PacBio_samples_file);
if($corrected_reads_file){$corrected_reads_file = abs_path($corrected_reads_file);}
if($RNAseq_samples_file){$RNAseq_samples_file = abs_path($RNAseq_samples_file);}
if($group_file){$group_file = abs_path($group_file);}
if($group_compare_file){$group_compare_file = abs_path($group_compare_file);}
$parameters_file = abs_path($parameters_file);
$database_file = abs_path($database_file);

if (!$project_name || !$species || !$category || !$PacBio_samples_file || !$parameters_file || !$database_file|| !$prefix || !$output_dir || $help){    
	die "$usage\n";
}
die "start-step must bigger than the end-step\n" if ($step >$cut);
#my @skips=("False","False","False","False","False","False","False","False","False","False","False","False");
my @skips=(0,0,0,0,0,0,0,0,0,0,0,0);
my @excls=split /,/,$skip;
foreach my $excl (@excls){$skips[$excl-1]=1;}
#####################################################################################
#                       write project infomation                                    #
##################################################################################### 
system("mkdir $output_dir") unless (-e "$output_dir");
system("mkdir $output_dir/result") unless (-e "$output_dir/result");
open  PI,">$output_dir/result/base_info.json" or die;
print PI "{\"projectName\":\"$project_name\" , \"species\":\"$species\"}\n";
close PI;
#####################################################################################
#                     Initialize the parameters                                     #
##################################################################################### 
open PARAMETERS, "$parameters_file" or die "Cannot open parameters_file: $!";
my %parameters;
while (<PARAMETERS>) {
	chomp;
	s/\s+$//;
	s/\r$//;
	next if (/^\#/ or /^\s*$/);
	my @arr = split /\s+/;
	$parameters{$arr[0]} = $arr[1];
}
close PARAMETERS;
#parameters for IsoSeq
my $minfull = $parameters{"minfull"} // 0;
my $minpredit = $parameters{"minpredit"} // 75;
my $minseq = $parameters{"minseq"} //300;
my $hq = $parameters{"hq"} //0.99;
#parameters for LSC
#parameters for CD-hit-est
my $c=$parameters{"c"} //0.9;
my $T=$parameters{"T"} //1;
my $G=$parameters{"G"} //1;
my $U=$parameters{"U"} //99999999;
my $s=$parameters{"s"} //0.0;
my $d=$parameters{"d"} //20;
my $p=$parameters{"p"} //0;
##parameters for Annotation
my $blast_cpu = $parameters{"blast_cpu"} // 50;
my $blast_e = $parameters{"blast_e"} // 1e-5;
my $blast_p = $parameters{"blast_p"} // "blastx";
my $blast_cut = $parameters{"blast_cut"} // 200;
my @databases;
if ($parameters{"cog"}) {push @databases, "--cog";}
if ($parameters{"kog"}) {push @databases, "--kog";}
if ($parameters{"nr"}) {push @databases, "--nr";}
if ($parameters{"nt"}) {push @databases, "--nt";}
if ($parameters{"swissprot"}) {push @databases, "--swissprot";}
if ($parameters{"tremble"}) {push @databases, "--tremble";}
if ($parameters{"kegg"}) {push @databases, "--kegg";}
if ($parameters{"GO"}) {push @databases, "--GO";}
if ($parameters{"all"}) {push @databases, "--all";}
my $data_base_choose = join " ", @databases;
#parameters for Differential expression
my $P = $parameters{"P"} // 0.05;
my $C = $parameters{"C"} // 2;
my $fdr = $parameters{"fdr"} // 0.05;
my $foldchange = $parameters{"foldchange"} // 1;

#####################################################################################
#              Step 1, obtain full length transcripts using IsoSeq                  #
##################################################################################### 
	system("mkdir $output_dir/worksh") unless (-e "$output_dir/worksh");
	system("perl $Bin/IsoSeq/H5_input.pl -infile $PacBio_samples_file");
	open SAMPLES, "$PacBio_samples_file.2" or die "Cannot open PacBio_samples_file: $!";
while (<SAMPLES>) {
	chomp;
	my @array = split /\t/, $_;
	open  SH,">$output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/IsoSeq.sh";
	print SH "echo ===start at : `date` ===\n";
	if($Iso_method eq "Web"){
	print SH "perl $Bin/IsoSeq/smrtanalysis_set.pl --input $Bin/IsoSeq/settings.xml --minfull $minfull --minpredit $minpredit --minseq $minseq --cDNA $array[1] --hq $hq --output settings.xml\n$Bin/01assemble/fofnToSmrtpipeInput.py input.fofn >input.xml\n. \$(/share/public/software/smrtanalysis_2.3.0/install/smrtanalysis_2.3.0.140936/admin/bin/getsetupfile) && /share/public/software/smrtanalysis_2.3.0/install/smrtanalysis_2.3.0.140936/admin/bin/checktmpdir --verbose --create --dir $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/tmp && export MPLCONFIGDIR=/share/public/software/smrtanalysis_2.3.0/userdata/jobs/.matplotlib && ([[ ! -e /share/public/software/smrtanalysis_2.3.0/userdata/jobs/.matplotlib ]] && mkdir /share/public/software/smrtanalysis_2.3.0/userdata/jobs/.matplotlib ) || true && smrtpipe.py -D TMP=$output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/tmp --nohtml  --distribute  --params=settings.xml --output=$output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1] xml:input.xml 2>smrtpipe.stderr 1> smrtpipe.stdout && \\\n";}
	else{
	print SH "if [ ! -d $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results ]; then\n mkdir $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results\n fi && \\\n";
	print SH "source /share/public/software/smrtanalysis_2.3.0/install/smrtanalysis_2.3.0.140936/etc/setup.sh && \\\n";
	print SH "ConsensusTools.sh CircularConsensus --minFullPasses 0 --minPredictedAccuracy 75 -p /share/public/software/smrtanalysis_2.3.0/install/smrtanalysis_2.3.0.140936/analysis/etc/algorithm_parameters/2015-11 --numThreads 30 --fofn $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/input.fofn -o $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data && cat $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/*.ccs.fasta > $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/reads_of_insert.fasta && cat $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/*.ccs.fastq > $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/reads_of_insert.fastq && ls $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/*.ccs.h5 > $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/reads_of_insert.fofn && reads_of_insert_report.py --output-dir $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/reads_of_insert.fofn $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results/reads_of_insert_report.json && \\\n";
	print SH "pbtranscript.py classify $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/reads_of_insert.fasta $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/isoseq_draft.fasta --cpus 30 --min_seq_len 300 --flnc $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/isoseq_flnc.fasta --nfl $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/isoseq_nfl.fasta && isoseq_classify_report.py --output-dir $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/isoseq_flnc.fasta $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/isoseq_draft.classify_summary.txt $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results/isoseq_classify.json && \\\n";
	print SH "pbtranscript.py cluster $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/isoseq_flnc.fasta $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/consensus_isoforms.fasta --nfl_fa $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/isoseq_nfl.fasta -d $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data --ccs_fofn $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/reads_of_insert.fofn --bas_fofn $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/input.fofn --quiver --blasr_nproc 30 --quiver_nproc 20 --cDNA_size $array[1] --hq_quiver_min_accuracy 0.99 --hq_isoforms_fa $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/polished_high_qv_consensus_isoforms.fasta --hq_isoforms_fq $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/polished_high_qv_consensus_isoforms.fastq --lq_isoforms_fa $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/polished_low_qv_consensus_isoforms.fasta --lq_isoforms_fq $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/polished_low_qv_consensus_isoforms.fastq && isoseq_cluster_report.py --output-dir $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/consensus_isoforms.fasta $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/data/cluster_summary.txt $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/results/isoseq_cluster.json && \\\n";
	}
#	print SH "source /home/xingsk/.bashrc && \\\n";
	print SH "python $Bin/IsoSeq/iso_seq_id.py ./data/ $array[0] $array[1] ./data/ && \\\n";#add sample alias to transcript name
	print SH "perl $Bin/IsoSeq/classify_pie.pl -in ./data/isoseq_draft.classify_summary.txt -out ./results/classify_summary && \\\n";
	print SH "python $Bin/IsoSeq/stat_length.py --in ./data/reads_of_insert.fasta --out ./results/roi_length.xls --cDNA $array[1] && \\\n";
	print SH "python $Bin/IsoSeq/stat_length.py --in ./data/isoseq_flnc.fasta --out ./results/flnc_length.xls --cDNA $array[1] && \\\n";
	print SH "python $Bin/IsoSeq/stat_length.py --in ./data/consensus_isoforms.fasta --out ./results/consensus_length.xls --cDNA $array[1] && \\\n";
	print SH "cp PacBio_rawdata*png ./results && \\\n";
        print SH "echo \"IsoSeq finish\" >IsoSeq.mark && \\\n";
	print SH "cp IsoSeq.mark $output_dir/worksh/$array[1]_$array[0].IsoSeq.mark || (touch $output_dir/worksh/IsoSeq.$array[0]_$array[1].failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/IsoSeq/$array[0]/$array[0]_$array[1]/IsoSeq.sh $output_dir/worksh/IsoSeq.$array[0]_$array[1].sh");
if($step<2 && $cut>0){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir -sample $array[0]/$array[0]_$array[1] -step 1 -cut 1 1>$output_dir/worksh/qsub.IsoSeq.$array[0]_$array[1].sh.o 2>$output_dir/worksh/qsub.IsoSeq.$array[0]_$array[1].sh.e &");
	}
}
	close SAMPLES;
#####################################################################################
#      Step 2,  full length transctipts correction with RNA-seq data using LSC      #
##################################################################################### 
	open SAMPLES, "$PacBio_samples_file.3" or die "Cannot open PacBio_samples_file: $!";
	my $j=0;
while (<SAMPLES>) {
	$j++;
	chomp;
	my @array=split /\t/,$_;
	unless($skips[1]){
if(-e "$corrected_reads_file") {
	system("mkdir $output_dir/result/LSC") unless (-e "$output_dir/result/LSC");
	system("mkdir $output_dir/result/LSC/$array[0]") unless (-e "$output_dir/result/LSC/$array[0]");
	open  (SH, "$corrected_reads_file");
	my %hash=();
	my $key;
	while (my $line=<SH>)
	{
	chomp $line;
	my @lsc=split /\t/,$line;
	$key=$lsc[0];
	shift @lsc;
	$hash{$key}=join(" --short_reads ",@lsc);
	}
	close SH;
	if($hash{$array[0]}){
	open  SH,">$output_dir/result/LSC/$array[0]/LSC.sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "cat $output_dir/result/IsoSeq/$array[0]/*/data/low_qv_consensus_isoforms.fasta >cat_lq_consensus_isoforms.fasta && \\\n";
	print SH "python /share/public/software/LSC/LSC-2.0/bin/runLSC.py --long_reads cat_lq_consensus_isoforms.fasta --short_reads $hash{$array[0]} -o $output_dir/result/LSC/$array[0] --short_read_file_type fq --tempdir $output_dir/result/LSC/$array[0]/tmp --threads 15 && \\\n";
	print SH "echo \"LSC finish\" >LSC.mark && \\\n";
	print SH "cp LSC.mark $output_dir/worksh/LSC.$array[0].mark || (touch $output_dir/worksh/LSC.$array[0].failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/LSC/$array[0]/LSC.sh $output_dir/worksh/LSC.$array[0].sh");
if($step<3 && $cut>1){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir -sample $array[0] -fileno $array[1] -step 2 -cut 2 1>$output_dir/worksh/qsub.LSC.$array[0].sh.o 2>$output_dir/worksh/qsub.LSC.$array[0].sh.e &");
		}
	}else{system("rm -rf $output_dir/result/LSC/$array[0]");}	
	}
}
#####################################################################################
#          Step 3,  remove redundant full length transcripts using CD-hit-est       #
##################################################################################### 
	system("mkdir $output_dir/result/CD-hit-est") unless (-e "$output_dir/result/CD-hit-est");
	system("mkdir $output_dir/result/CD-hit-est/$array[0]") unless (-e "$output_dir/result/CD-hit-est/$array[0]");
	open  SH,">$output_dir/result/CD-hit-est/$array[0]/CD-hit-est.sh";
	print SH "echo ===start at : `date` ===\n";
	unless (-e "$output_dir/result/LSC/$array[0]"){
	print SH "cat $output_dir/result/IsoSeq/$array[0]/*/data/high_qv_consensus_isoforms.fasta >cat_consensus_isoforms.fasta && \\\n";
	print SH "/share/public/software/CD-HIT/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i cat_consensus_isoforms.fasta -o UniIso.fa -c $c -T $T -G $G -U $U -s $s -d $d -p $p && \\\n";}
	else {
	print SH "cat $output_dir/result/IsoSeq/$array[0]/*/data/high_qv_consensus_isoforms.fasta >cat_hq_consensus_isoforms.fasta && \\\n";
	print SH "cat cat_hq_consensus_isoforms.fasta $output_dir/result/LSC/$array[0]/full_LR.fa >cat_consensus_isoforms.fasta && \\\n";
	print SH "/share/public/software/CD-HIT/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i cat_consensus_isoforms.fasta  -o UniIso.fa -c $c -T $T -G $G -U $U -s $s -d $d -p $p && \\\n";}
	print SH "echo \"CD-hit-est finish\" >CD-hit-est.mark && \\\n";
	print SH "cp CD-hit-est.mark $output_dir/worksh/$array[0].CD-hit-est.mark || (touch $output_dir/worksh/CD-hit-est.$array[0].failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/CD-hit-est/$array[0]/CD-hit-est.sh $output_dir/worksh/CD-hit-est.$array[0].sh");
if ($step<4 && $cut>2){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir  -sample $array[0] -fileno $array[1] -step 3  -cut 3 1>$output_dir/worksh/qsub.CD-hit-est.$array[0].sh.o 2>$output_dir/worksh/qsub.CD-hit-est.$array[0].sh.e &");
	}	
#####################################################################################
#                       Step 4, Functional  annotation                              #
##################################################################################### 
	unless($skips[3]){
	system("mkdir $output_dir/result/Annotation") unless (-e "$output_dir/result/Annotation");
	system("mkdir $output_dir/result/Annotation/$array[0]") unless (-e "$output_dir/result/Annotation/$array[0]");
	open  SH,">$output_dir/result/Annotation/$array[0]/Annotation.sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "perl $Bin/annotation/Gene_Func_Anno_Pipline.pl --path $output_dir/result/CD-hit-est/$array[0]/ $data_base_choose --Blastp $blast_p --Blast_cpu $blast_cpu --Blast_e $blast_e --Blast_cut $blast_cut --Database_cfg $database_file --od $output_dir/result/Annotation/$array[0]/ --species $category && \\\n";
	print SH "echo \"Annotation finish\" >Annotation.mark && \\\n";
	print SH "cp Annotation.mark $output_dir/worksh/Annotation.$array[0].mark || (touch $output_dir/worksh/Annotation.$array[0].failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/Annotation/$array[0]/Annotation.sh $output_dir/worksh/Annotation.$array[0].sh");
if ($step<5 && $cut>3){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir  -sample $array[0]  -step 4  -cut 4 1>$output_dir/worksh/qsub.Annotation.$array[0].sh.o 2>$output_dir/worksh/qsub.Annotation.$array[0].sh.e &");
	}
}
#####################################################################################
#                          Step 5, CDS prediction using Transdecoder                #
##################################################################################### 
	unless($skips[4]){
	system("mkdir $output_dir/result/Transdecoder") unless (-e "$output_dir/result/Transdecoder");
	system("mkdir $output_dir/result/Transdecoder/$array[0]") unless (-e "$output_dir/result/Transdecoder/$array[0]");
	open  SH,">$output_dir/result/Transdecoder/$array[0]/Transdecoder.sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "ln -s $output_dir/result/CD-hit-est/$array[0]/UniIso.fa UniIso.fa && \\\n";
	print SH "/share/public/software/TransDecoder/TransDecoder-3.0.0/TransDecoder.LongOrfs -t UniIso.fa && \\\n";
	print SH "/share/public/software/ncbi-blast-2.2.25+/bin/blastp -query ./UniIso.fa.transdecoder_dir/longest_orfs.pep -db /share/public/database/Uniprot/2016-08-03/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 >blastp.outfmt6 && \\\n";
	print SH "/share/public/software/hmmer/hmmer-3.1b2-linux-intel-x86_64/bin/hmmscan --cpu 8 --domtblout pfam.domtblout /share/public/software/pfam/test/Pfam-A.hmm ./UniIso.fa.transdecoder_dir/longest_orfs.pep && \\\n";
	print SH "/share/public/software/TransDecoder/TransDecoder-3.0.0/TransDecoder.Predict -t UniIso.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_orf && \\\n";
	print SH "perl $Bin/Transdecoder_sort.pl UniIso.fa.transdecoder.cds && \\\n";
	print SH "perl $Bin/Transdecoder_sort.pl UniIso.fa.transdecoder.pep && \\\n";
	print SH "echo \"Transdecoder finish\" >Transdecoder.mark && \\\n";
	print SH "cp Transdecoder.mark $output_dir/worksh/Transdecoder.$array[0].mark || (touch $output_dir/worksh/Transdecoder.$array[0].failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/Transdecoder/$array[0]/Transdecoder.sh $output_dir/worksh/Transdecoder.$array[0].sh");
if ($step<6 && $cut>4){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir  -sample $array[0]  -step 5  -cut 5 1>$output_dir/worksh/qsub.Transdecoder.$array[0].sh.o 2>$output_dir/worksh/qsub.Transdecoder.$array[0].sh.e &");
	}
}
#####################################################################################
#                   Step 6, Gene family identification using OrthoMCL               #
##################################################################################### 
	unless($skips[5]){
	system("mkdir $output_dir/result/OrthoMCL") unless (-e "$output_dir/result/OrthoMCL");
	system("mkdir $output_dir/result/OrthoMCL/$array[0]") unless (-e "$output_dir/result/OrthoMCL/$array[0]");
	open  SH,">$output_dir/result/OrthoMCL/$array[0]/Protein.config";
	print SH "$array[0]	$output_dir/result/Transdecoder/$array[0]/UniIso.fa.transdecoder.pep	7\n";
	close SH;
	open  SH,">$output_dir/result/OrthoMCL/$array[0]/OrthoMCL.sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "perl $Bin/orthomcl.main_pipe.v2.pl -cfg Protein.config -o $output_dir/result/OrthoMCL/$array[0]/ -prefix $prefix && \\\n";
	print SH "make -j 10 -f orthomcl.make && \\\n";
	print SH "echo \"OrthoMCL finish\" >OrthoMCL.mark && \\\n";
	print SH "cp OrthoMCL.mark $output_dir/worksh/OrthoMCL.$array[0].mark && cp OrthoMCL.mark $output_dir/worksh/OrthoMCL.mark.$j || (touch $output_dir/worksh/OrthoMCL.$array[0].failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/OrthoMCL/$array[0]/OrthoMCL.sh $output_dir/worksh/OrthoMCL.$array[0].sh");
if ($step<7 && $cut>5){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir -sample $array[0] -fileno $j -step 6  -cut 6 1>$output_dir/worksh/qsub.OrthoMCL.$array[0].sh.o 2>$output_dir/worksh/qsub.OrthoMCL.$array[0].sh.e &");
	}
}
#####################################################################################
#                                Step 7, SSR prediction                             #
##################################################################################### 
	unless($skips[6]){
	system("mkdir $output_dir/result/SSR") unless (-e "$output_dir/result/SSR");
	system("mkdir $output_dir/result/SSR/$array[0]") unless (-e "$output_dir/result/SSR/$array[0]");
        open  SH,">$output_dir/result/SSR/$array[0]/SSR.sh";
        print SH "echo ===start at : `date` ===\n";
#        print SH "perl $Bin/SSR/SSR_analysis.pl -fasta $output_dir/result/CD-hit-est/$array[0]/UniIso.fa -outdir $output_dir/result/SSR/$array[0] -name $prefix && \\\n";
        print SH "$Bin/SSR/runSSR_primer -a $output_dir/result/CD-hit-est/$array[0]/UniIso.fa -dir $output_dir/result/SSR/$array[0] -name $prefix && \\\n";
        print SH "sh runSSR_primer.sh && \\\n";
        print SH "echo \"SSR finish\" >SSR.mark && \\\n";
	print SH "cp SSR.mark $output_dir/worksh/SSR.$array[0].mark || (touch $output_dir/worksh/SSR.$array[0].failed;touch $output_dir/worksh/failed) \n";
        print SH "echo ===end at : `date` ===\n";
        close SH;
system("cp $output_dir/result/SSR/$array[0]/SSR.sh $output_dir/worksh/SSR.$array[0].sh");
if ($step<8 && $cut>6){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir  -sample $array[0]  -step 7  -cut 7 1>$output_dir/worksh/qsub.SSR.$array[0].sh.o 2>$output_dir/worksh/qsub.SSR.$array[0].sh.e &");
        }
}
#####################################################################################
#                       Step 8, LncRNA prediction                                   #
##################################################################################### 
	unless($skips[7]){
	system("mkdir $output_dir/result/LncRNA") unless (-e "$output_dir/result/LncRNA");
	system("mkdir $output_dir/result/LncRNA/$array[0]") unless (-e "$output_dir/result/LncRNA/$array[0]");
        open  SH,">$output_dir/result/LncRNA/$array[0]/LncRNA.sh";
        print SH "echo ===start at : `date` ===\n";
	print SH "source /share/public/pipeline/RNA/config/bashrc_RNA_seq_pacbio \n";
        print SH "ln -s $output_dir/result/CD-hit-est/$array[0]/UniIso.fa UniIso.fa && \\\n";
        print SH "/share/public/software/LncRNAs/lncrnas-pipeline/ncrna_pipeline -f UniIso.fa -p 10 -c plek && \\\n";
        print SH "python  $Bin/IsoSeq/stat_length.py --in UniIso_ncrna.fa --out lnc_length.xls --cDNA cDNA && \\\n";
        print SH "source ~/.bashrc \n";
	print SH "Rscript $Bin/lnc/lnc_length.R lnc_length.xls lncRNA_length_cout && \\\n";
        print SH "echo \"LncRNA finish\" >LncRNA.mark && \\\n";
	print SH "cp LncRNA.mark $output_dir/worksh/LncRNA.$array[0].mark || (touch $output_dir/worksh/LncRNA.$array[0].failed;touch $output_dir/worksh/failed) \n";
        print SH "echo ===end at : `date` ===\n";
        close SH;
system("cp $output_dir/result/LncRNA/$array[0]/LncRNA.sh $output_dir/worksh/LncRNA.$array[0].sh");
if ($step<9 && $cut>7){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir  -sample $array[0]  -step 8  -cut 8 1>$output_dir/worksh/qsub.LncRNA.$array[0].sh.o 2>$output_dir/worksh/qsub.LncRNA.$array[0].sh.e &");
        }
}
}
	close SAMPLES;
#####################################################################################
#                       Step 9,  Get reference isoform                              #
##################################################################################### 
unless($skips[8]){
if($RNAseq_samples_file){
	system("mkdir $output_dir/result/RefIso") unless (-e "$output_dir/result/RefIso");
	open  SH,">$output_dir/result/RefIso/RefIso.sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "until test \$(ls -l $output_dir/worksh/ | grep -c '.CD-hit-est.mark\$') = $j\ndo\n sleep 1m\ndone && \\\n";#match the *.CD-hit-est.mark files in the directory
	if($j==1){
        print SH "cat $output_dir/result/CD-hit-est/*/UniIso.fa >$output_dir/result/RefIso/RefIso.fa && \\\n";}
	else{
        print SH "cat $output_dir/result/CD-hit-est/*/UniIso.fa >$output_dir/result/RefIso/Multiple_isoforms.fasta && \\\n";
	print SH "/share/public/software/CD-HIT/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i Multiple_isoforms.fasta -o RefIso.fa -c $c -T $T -G $G -U $U -s $s -d $d -p $p && \\\n";}
	print SH "echo \"RefIso finish\" >RefIso.mark && \\\n";
	print SH "cp RefIso.mark $output_dir/worksh/RefIso.mark || (touch $output_dir/worksh/RefIso.failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/RefIso/RefIso.sh $output_dir/worksh/RefIso.sh");
if ($step<10 && $cut>8){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir  -sample $project_name -fileno $j -step 9  -cut 9 1>$output_dir/worksh/qsub.RefIso.sh.o 2>$output_dir/worksh/qsub.RefIso.sh.e &");
		}
	}
}	
#####################################################################################
#                       Step 10, Expression evaluation                              #
##################################################################################### 
my @RNAseq_samples;
my @RNAseq_sample_name;
my $RNAseq_sample_names;
	my $i=0;
unless($skips[9]){
if($RNAseq_samples_file) {
	open SAMPLES, "$RNAseq_samples_file" or die "Cannot open RNAseq_samples_file: $!";
	system("mkdir $output_dir/result/RSEM") unless (-e "$output_dir/result/RSEM");
	open SH,">$output_dir/result/RSEM/rsem-prepare-reference.sh";
	print SH "echo \"start at : `date` \" >reference_prepare_start.mark\n";
	print SH "echo ===start at : `date` ===\n";
	print SH "/share/public/software/trinity/trinityrnaseq_r20140717/trinity-plugins/rsem/rsem-prepare-reference --bowtie2 --no-polyA $output_dir/result/RefIso/RefIso.fa $prefix && \\\n";
	print SH "echo \"rsem-prepare-reference finish\" >reference_prepare.mark && \\\n";
	print SH "cp reference_prepare.mark $output_dir/worksh/reference_prepare.mark || (touch $output_dir/worksh/reference_prepare.failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/RSEM/rsem-prepare-reference.sh $output_dir/worksh/rsem-prepare-reference.sh");
while (<SAMPLES>) {
	chomp;
	s/\s+$//;
	s/\r$//;
	next if (/^\#/ or /^$/);
	@RNAseq_samples=split /\t/,$_;
	push @RNAseq_sample_name,$RNAseq_samples[0];
	open  SH,">$output_dir/result/RSEM/expression.$RNAseq_samples[0].sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "mkdir $RNAseq_samples[0] && \\\n";
	if($RNAseq_samples[2] eq "NA"){
	print SH "rsem-calculate-expression -p 8 --bowtie2 --estimate-rspd $RNAseq_samples[1]  $prefix $RNAseq_samples[0]/$RNAseq_samples[0] && \\\n";
	print SH "$Bin/expression/bam2map_stat $RNAseq_samples[0] ./$RNAseq_samples[0]/$RNAseq_samples[0].transcript.bam se ./$RNAseq_samples[0]/$RNAseq_samples[0].mapped.stat.xls && \\\n";
	print SH "awk '{print \$1 \"\\t\" \$3 \"\\t\" \$5 \"\\t\" \$7}' ./$RNAseq_samples[0]/$RNAseq_samples[0].isoforms.results > ./$RNAseq_samples[0]/$RNAseq_samples[0]_TranscriptExp.xls && \\\n";
	print SH "$Bin/expression/RSEM_result --sample $RNAseq_samples[0] --indir ./$RNAseq_samples[0] --outdir ./$RNAseq_samples[0] && \\\n";
	print SH "cat ./$RNAseq_samples[0]/$RNAseq_samples[0].mapped.stat.xls | sed -n '2p' > ./$RNAseq_samples[0]/$RNAseq_samples[0].mapped.stat.2.xls && \\\n";}
	else{	
	print SH "rsem-calculate-expression -p 8 --paired-end --bowtie2 --estimate-rspd $RNAseq_samples[1]  $RNAseq_samples[2] $prefix $RNAseq_samples[0]/$RNAseq_samples[0] && \\\n";
	print SH "$Bin/bam2map_stat $RNAseq_samples[0] ./$RNAseq_samples[0]/$RNAseq_samples[0].transcript.bam pe ./$RNAseq_samples[0]/$RNAseq_samples[0].mapped.stat.xls && \\\n";
	print SH "awk '{print \$1 \"\\t\" \$3 \"\\t\" \$5 \"\\t\" \$7}' ./$RNAseq_samples[0]/$RNAseq_samples[0].isoforms.results > ./$RNAseq_samples[0]/$RNAseq_samples[0]_TranscriptExp.xls && \\\n";
	print SH "$Bin/expression/RSEM_result --sample $RNAseq_samples[0] --indir ./$RNAseq_samples[0] --outdir ./$RNAseq_samples[0] && \\\n";
	print SH "cat ./$RNAseq_samples[0]/$RNAseq_samples[0].mapped.stat.xls | sed -n '2p' > ./$RNAseq_samples[0]/$RNAseq_samples[0].mapped.stat.2.xls && \\\n";}
	print SH "echo \"evaluate expression finish\" >$RNAseq_samples[0].expression.mark && \\\n";
	print SH "cp $RNAseq_samples[0].expression.mark $output_dir/worksh/expression.$RNAseq_samples[0].mark || (touch $output_dir/worksh/expression.$RNAseq_samples[0].failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/RSEM/expression.$RNAseq_samples[0].sh $output_dir/worksh/expression.$RNAseq_samples[0].sh");
if ($step<11 && $cut>9){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir  -sample $RNAseq_samples[0]  -step 10  -cut 10 1>$output_dir/worksh/qsub.expression.$RNAseq_samples[0].sh.o 2>$output_dir/worksh/qsub.expression.$RNAseq_samples[0].sh.e &");
	}
	$i++; 
	sleep 20;#wait for the file reference_prepare_start.mark to be produced by rsem-prepare-reference.sh of first Illumina sample   
}	
	close SAMPLES;
} 
}
	$RNAseq_sample_names=join ",",@RNAseq_sample_name;
	open SH,">$output_dir/result/RSEM/sample_list.txt";
	print SH $RNAseq_sample_names."\n";
	close SH; 
#####################################################################################
#                       Step 11, Differential expression                            #
##################################################################################### 
unless($skips[10]){
if($group_compare_file){
	system("mkdir $output_dir/result/edgeR") unless (-e "$output_dir/result/edgeR");
	open  SH,">$output_dir/result/edgeR/diff_exp.sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "cat $output_dir/result/Annotation/*/Result/UniIso.fa.GO.list.xls > RefIso.fa.GO.list.xls && \\\n";
	print SH "cat $output_dir/result/Annotation/*/Result/UniIso.fa.KEGG.blast.xls > RefIso.fa.KEGG.blast.xls && \\\n";
	print SH "perl $Bin/differential_expression/abundance_estimats_to_matrix.pl --est_method RSEM --name_sample $RNAseq_sample_names --out_prefix $prefix.genes --exp_sample $output_dir/result/RSEM && \\\n";
	print SH "perl $Bin/differential_expression/run_DE_analysis.pl --matrix $prefix.genes.counts.matrix --method $DE_method --samples_file $group_file --contrasts $group_compare_file --FDR $fdr --log2FC $foldchange && \\\n";
	print SH "perl $Bin/differential_expression/analyze_diff_expr.pl --matrix $output_dir/result/edgeR/$prefix.genes.TMM.fpkm.matrix -P $P -C $C --samples $group_file && \\\n";
#	print SH "perl $Bin/differential_expression/diff_exp_table_GO_kegg.pl -matrix $output_dir/result/edgeR/$prefix.genes.TMM.fpkm.matrix -dir $output_dir/result/edgeR -sample $prefix -fdr $fdr -foldchange $foldchange && \\\n";
	print SH "export PYTHONPATH=/share/public/software/kobas-3.0//src/:/share/public/software/Python-2.7.13/lib/python2.7/site-packages/:/share/public/software/Python-2.7.13/lib/python2.7/site-packages/:/share/public/software/WDLIB/lib/FALCON_unzip_0619.2017/lib/python2.7/site-packages/:/share/public/software/falcon_unzip/lib/python2.7/site-packages/:/share/public/software/kobas-3.0/src/:/share/public/software/falcon_unzip/lib/python2.7/site-packages/ \n";
	print SH "perl $Bin/differential_expression/diff_exp_table_repeat_GO_kegg.pl -repeat $group_file -matrix $output_dir/result/edgeR/$prefix.genes.TMM.fpkm.matrix -dir $output_dir/result/edgeR -sample $prefix -fdr $fdr -foldchange $foldchange && \\\n";
	print SH "echo \"differential expression finish\" >diff_exp.mark && \\\n";
	print SH "cp diff_exp.mark $output_dir/worksh/diff_exp.mark || (touch $output_dir/worksh/diff_exp.failed;touch $output_dir/worksh/failed) \n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
system("cp $output_dir/result/edgeR/diff_exp.sh $output_dir/worksh/diff_exp.sh");
if ($step<12 && $cut>10){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir -sample $project_name -fileno $i -step 11  -cut 11 1>$output_dir/worksh/qsub.diff_exp.sh.o 2>$output_dir/worksh/qsub.diff_exp.sh.e &");
	}
}
}
#####################################################################################
#                     Step 12, Package result and web report                        #
##################################################################################### 
	open  SH,">$output_dir/Web_report.sh";
	print SH "echo ===start at : `date` ===\n";
	print SH "perl $Bin/package_result.pl $output_dir $prefix $PacBio_samples_file.3 || touch Web_report.failed \n";
	print SH "perl $Bin/gene_anno.pl -indir $output_dir/PacBio_RNA_denovo_report/result/advanced/Differential_expression -ann $output_dir/PacBio_RNA_denovo_report/result/primary/Annotation/*/Integrated_Function.annotation.xls";
	print SH "echo \"all finished\" >all.finish.mark\n";
	print SH "echo ===end at : `date` ===\n";
	close SH;
if ($step<13 && $cut>11){
system ("nohup perl $Bin/qsub.pl -inputdir $output_dir -sample $project_name  -step 12  -cut 12 1>qsub.Web_report.sh.o 2>qsub.Web_report.sh.e &"); }
