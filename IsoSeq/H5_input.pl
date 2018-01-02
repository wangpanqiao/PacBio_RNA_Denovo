# Writer:       Xing Shilai 
# # Program Date:   2016.08.15
#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Getopt::Long;
my $usage =
"
Usage:

  Options:

  -h                     Help
  -infile       <str>    the full path of the file consisting of six coloums,<Sample name>,<cDNA>,<*p0.1.bax.h5>,<*p0.2.bax.h5>,<*p0.3.bax.h5>,<Alias> 
";

my($infile, $sample, $cDNA, $out, $content, $help);
GetOptions(
          "infile=s" => \$infile,
          "help|h=s" => \$help
);
if(!$infile || $help){
        die "$usage\n";
}
$out=dirname($infile);
system("mkdir $out/result/IsoSeq") unless (-e "$out/result/IsoSeq");
open SAMPLES, "$infile" or die "Cannot open PacBio_samples_file: $!";
while (<SAMPLES>) {
        chomp;
        s/\s+$//;
        s/\r$//;
        next if (/^\#/ or /^$/);
        my @array = split /\t/, $_;
        if(-e "$out/result/IsoSeq/$array[0]/$array[0]_$array[1]/input.fofn") {system("rm -r $out/result/IsoSeq/$array[0]/$array[0]_$array[1]/input.fofn");}
}
close SAMPLES;
open SAMPLES, "$infile";
open OUT, ">PacBio_sample.xls";
print OUT "样品名"."\t"."Cell编号"."\t"."文库类型"."\t"."Polymerase Read Bases"."\t"."Polymerase Reads"."\t"."Polymerase Read N50"."\t"."Polymerase Read Length"."\t"."Polymerase Read Quality"."\n";
my $line;
my @part1=();
my @part2=();
my $i;
while (<SAMPLES>) {
        $i++;
	chomp;
        s/\s+$//;
        s/\r$//;
        next if (/^\#/ or /^$/);
        $line=$_;
        my @array = split /\t/, $line;
	my @cell = split /\//,$array[2];
	my @dot = split /\./,$cell[8];
	print OUT $array[0]."\t".$cell[6]."\t".$array[1]."\n";
        $line=$array[0]."\t".$array[1]."\t".$array[5];
        push @part1,$line;
        system("mkdir $out/result/IsoSeq/$array[0]") unless (-e "$out/result/IsoSeq/$array[0]");
        system("mkdir $out/result/IsoSeq/$array[0]/$array[0]_$array[1]") unless (-e "$out/result/IsoSeq/$array[0]/$array[0]_$array[1]");
        open PACBIO, ">>$out/result/IsoSeq/$array[0]/$array[0]_$array[1]/input.fofn" or die "Cannot open input.fofn: $!";
        print PACBIO "$array[2]\n$array[3]\n$array[4]\n";
        close PACBIO;
	if($array[2]=~/(.*)\/Analysis_Results\/.*/){
system("tree --charset X $1 -o $out/result/IsoSeq/$array[0]/$array[0]_$array[1]/PacBio_rawdata.$i.txt");
system("mogrify -format png -density 1024x1024 -trim $out/result/IsoSeq/$array[0]/$array[0]_$array[1]/PacBio_rawdata.$i.txt");
}
}
close SAMPLES;
close OUT;
my %count;
my @uniq = grep { ++$count{ $_ } < 2; } @part1;
        open PACBIO, ">$infile.2";
foreach $line(@uniq)
{ 
        print PACBIO "$line\n";
       my @array = split /\t/, $line;
	push @part2,$array[0];
}
        close PACBIO;
        open PACBIO, ">$infile.3";
my (%hash,$key,$value); foreach $line (@part2){$hash{$line}++;} while(($key, $value) = each(%hash)) {print PACBIO $key."\t".$value."\n";}
        close PACBIO;


