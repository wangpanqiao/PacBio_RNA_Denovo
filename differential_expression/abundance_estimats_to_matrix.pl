#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin $Script);
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Data::Dumper;

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --est_method <method>  sample1.results sample2.results ...
# Required:
#
#  --name_sample <string>            name sample compatible with exp_sample files order, separated by comma;
#  --exp_sample  <files_dir>         files of abundance estimation results (RSEM or eXpress), separated by comma;
#
#  --est_method <string>           RSEM|eXpress  (needs to know what format to expect)
#
# Options:
#
#  --cross_sample_fpkm_norm <string>    TMM|UpperQuartile|none   (default: TMM)
#
#  --out_prefix <string>                default: 'matrix'
#
############################################################


__EOUSAGE__

    ;


my $help_flag;
my $est_method;
my $cross_sample_fpkm_norm = "TMM";
my $name_sample;
my $out_prefix = "matrix";
my $exp_sample;

&GetOptions('help|h' => \$help_flag,
            'name_sample=s' => \$name_sample,
            'exp_sample=s' => \$exp_sample,

            'est_method=s' => \$est_method,
            
            'cross_sample_fpkm_norm=s' => \$cross_sample_fpkm_norm,
            'out_prefix=s' => \$out_prefix,
            
            
            );


unless ($est_method && $name_sample && $exp_sample) {
    die $usage;
}

unless ($est_method =~ /^(RSEM|eXpress)$/i) {
    die "Error, dont recognize --est_method $est_method ";
}
unless ($cross_sample_fpkm_norm =~ /^(TMM|UpperQuartile|none)$/i) {
    die "Error, dont recognize --cross_sample_fpkm_norm $cross_sample_fpkm_norm ";
}

my @filenames = split/,/,$name_sample;
my @genes_results = `ls $exp_sample/\*/\*genes.results`;
$exp_sample = &format(@genes_results);
my @files = split/,/,$exp_sample;
print "$exp_sample\n";
print Dumper @files;
for (0..$#files) {
	print "$files[$_]\n";
	$files[$_] = &ABSOLUTE_DIR ($files[$_]);
}


=data_formats

## RSEM:

0       transcript_id
1       gene_id
2       length
3       effective_length
4       expected_count
5       TPM
6       FPKM
7       IsoPct


## eXpress:

0       bundle_id
1       target_id
2       length
3       eff_length
4       tot_counts
5       uniq_counts
6       est_counts
7       eff_counts
8       ambig_distr_alpha
9       ambig_distr_beta
10      fpkm
11      fpkm_conf_low
12      fpkm_conf_high
13      solvable

=cut
    
    ;

my ($acc_field, $counts_field, $fpkm_field);

if ($est_method =~ /^rsem$/i) {
    $acc_field = 0;
    $counts_field = 4;
    $fpkm_field = 6;
}
elsif ($est_method =~ /^express$/i) {
    $acc_field = 1;
    $counts_field = 7;
    $fpkm_field = 10;
}
else {
    die "Error, dont recognize --est_method $est_method ";
}

main: {
    
    my %data;
    
    foreach my $file (@files) {
        print STDERR "-reading file: $file\n";
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>; # ignore it
        while (<$fh>) {
            chomp;
            
            my @x = split(/\t/);
            my $acc = $x[$acc_field];
            my $count = $x[$counts_field];
            my $fpkm = $x[$fpkm_field];
            
            $data{$acc}->{$file}->{count} = $count;
            $data{$acc}->{$file}->{fpkm} = $fpkm;
        }
        close $fh;
    }
    
    print STDERR "\n\n* Outputting combined matrix.\n\n";
    
    my $counts_matrix_file = "$out_prefix.counts.matrix";
    my $fpkm_matrix_file = "$out_prefix.not_cross_norm.fpkm.tmp";
    open (my $ofh_counts, ">$counts_matrix_file") or die "Error, cannot write file $counts_matrix_file";
    open (my $ofh_fpkm, ">$fpkm_matrix_file") or die "Error, cannot write file $fpkm_matrix_file";

    { # check to see if they're unique
        my %filename_map = map { + $_ => 1 } @filenames;
        if (scalar keys %filename_map != scalar @filenames) {
            die "Error, the column headings: @filenames are not unique. please check your $name_sample and $exp_sample parameters?";
        }
    }
    

    print $ofh_counts "#geneID";
    print $ofh_fpkm "#geneID";
    print $ofh_counts join("\t", "", @filenames) . "\n";
    print $ofh_fpkm join("\t", "", @filenames) . "\n";

    foreach my $acc (keys %data) {
        
        print $ofh_counts "$acc";
        print $ofh_fpkm "$acc";
        
        foreach my $file (@files) {

            my $count = $data{$acc}->{$file}->{count};
            unless (defined $count) {
                $count = "NA";
            }
            my $fpkm = $data{$acc}->{$file}->{fpkm};
            unless (defined $fpkm) {
                $fpkm = "NA";
            }

            print $ofh_counts "\t$count";
            print $ofh_fpkm "\t$fpkm";
        }
        
        print $ofh_counts "\n";
        print $ofh_fpkm "\n";

    }
    close $ofh_counts;
    close $ofh_fpkm;
    
    if ($cross_sample_fpkm_norm =~ /^TMM$/i) {
        my $cmd = "perl $Bin/run_TMM_scale_matrix.pl --matrix $fpkm_matrix_file > $out_prefix.$cross_sample_fpkm_norm.fpkm.matrix";
        &process_cmd($cmd);
    }
    elsif ($cross_sample_fpkm_norm =~ /^UpperQuartile$/) {
        my $cmd = "perl $Bin/run_UpperQuartileNormalization_matrix.pl --matrix $fpkm_matrix_file > $out_prefix.$cross_sample_fpkm_norm.fpkm.matrix";
        &process_cmd($cmd);
    }
    elsif ($cross_sample_fpkm_norm =~ /^none$/i) {
        print STDERR "-not performing cross-sample normalization.\n";
    }
    
    print STDERR "Done.\n\n";
    
    exit(0);
}
    
####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR $cmd;
    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

sub format{
	my @array=@_;
	my @result;
	my $out;
	for (@array) {
		chomp;
		push @result, $_;
	}
	$out = join ",", @result;
}
