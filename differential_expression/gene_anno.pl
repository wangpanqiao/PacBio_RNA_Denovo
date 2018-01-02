#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;


my ($indir,$ann,$SNP,$AS);
GetOptions(
	'indir:s' => \$indir,
	'anno:s' => \$ann,
	'SNP' => \$SNP,
	'AS' => \$AS,
);

&USAGE unless ( $indir and $ann );

#############################################
my @FPKM = glob ("$indir/Quantification/*/*.FPKM.xls");
#my @DEG = glob ("$indir/DiffAnalysis/*_edgeR/*_{all_list,Up,Down}.xls");
my @DEG = glob ("$indir/*diff_exp.xls");
my @SNP = glob ("$indir/SNP/*_SNP/*variant_function.xls");
my @AS = glob ("$indir/Alternative_Splice/*_AS/*_altsplice.xls");

my %gene2ann;

#################gene annotation##############

open ANN,$ann or die $!;
while(<ANN>){
	chomp;
	next if /^#/;
	my @arr = split /\t/;
	my $gene = $arr[0];
#	splice(@arr,0,1);
	my $upper_gene = uc($gene);
	$gene2ann{$gene} = join ("\t",$arr[2],$arr[3],$arr[11],$arr[16],$arr[21]);
	$gene2ann{$upper_gene} = join ("\t",$arr[2],$arr[3],$arr[11],$arr[16],$arr[21]);

}
close ANN;

################Sample.RPKM.xls################

foreach ( @FPKM ) {
    if ( `grep KEGG_Pathway $_` ) {###judge the $genelist  whether annotated or not  
        print "The $_ had been annotated!\n";
        next;
    }
	&express_anno($_,"$_.ann","1");
}

################*_Up/Down/all_list.xls#########

foreach ( @DEG ) {
    if ( `grep KEGG_Pathway $_` ) {###judge the $genelist  whether annotated or not  
        print "The $_ had been annotated!\n";
        next;
    }
	&express_anno($_,"$_.ann","1");
}

###############*_altsplice.xls#################

if ( $AS ) {
    foreach ( @AS ) {
	if ( `grep KEGG_Pathway $_` ) {###judge the $in  whether annotated or not 
            print "The $_ had been annotated !\n";
            next;
        }

            &express_anno($_,"$_.ann",$AS);
    }
}

##############*variant_function.xls###########

if ( $SNP ) {
    foreach ( @SNP ) {
        if ( `grep KEGG_Pathway $_` ) {###judge the $in  whether annotated or not 
            print "The $_ had been annotated !\n";
	    next;
        }

    	    &snp_anno($_,"$_.ann",$SNP);
   }
}

###########################SUB FUNCTIONS################

sub express_anno {
	my ($genelist,$outfile,$type)=@_;
	
	open GENE, $genelist or die $!;
	open OUT, ">$outfile" or die $!;
	chomp (my $head = <GENE>);
	if ($head=~/GO_Anno/){
		print "Input file has been annotated\n";
		last;
	}
#	print OUT "$head\tGeneSymbol\tGeneID\tDescription\tBiological Process\tCellular Component\tMolecular Function\tKO_ID(KO_name)\tKEGG_Pathway\n";
	print OUT "$head\tGO_Anno\tKEGG_entry\tNR_annotation\tNT_annotation\tSwissprot_annotation\n";
	while(<GENE>){
            chomp;
            my $gene=(split/\t/)[0];
            $gene=(split /\+/,$gene)[0] if $gene =~/\+/;
            if($gene2ann{$gene}){
                print OUT "$_\t$gene2ann{$gene}\n";
            }else{
                print OUT "$_\t-\t-\t-\t-\t-\t-\t-\t-\n";
            }
       }
       close GENE;
       close OUT;
	
       ###rename *.ann files###
       if ( $type ) {
           my $dir = dirname $outfile;	
           `rename '.ann' '' $dir/*.ann`;
       }
}

sub uniq{
        my ($uniq)=@_;
        my %count;
        my @uniq=grep { ++$count{$_}< 2 } @{$uniq};
        return (@uniq);
}

sub snp_parse {
        my ($snp_type,$gene_name)=@_;
        my @gene_list=();

        if ( $snp_type =~ /intergenic/ ) {
             @gene_list=split (/,/,$gene_name);
             map { $_=~s/\(.*?\)// } @gene_list;

        }elsif ( $snp_type =~ /upstream|downstream|exonic|intronic|UTR3|UTR5/) {
             @gene_list=split (/,/,$gene_name);

        }elsif ( $snp_type =~ /splicing/ ) {
             @gene_list=split (/\),/,$gene_name);
             map { $_=~s/\(.*// } @gene_list;
        }
        return (@gene_list);
}


sub snp_anno{
	my ($in,$outfile,$type)=@_;
	
	if ( `grep KEGG_Pathway $in` ) {###judge the $in  whether annotated or not 
		print "The $in had been annotated !\n";
	}
		
	my $filename=basename $in;
	if( $filename =~ /\.variant_function\.xls/){
            open IN,$in or die $!;
	    open SNP_OUT,">$outfile" or die $!;
            chomp (my $head = <IN>);
            print SNP_OUT "$head\tGeneName(GeneSymbol)\tGeneID\tDescription\tBiological Process\tCellular Component\tMolecular Function\tKO_ID(KO_name)\tKEGG_Pathway\n";
            while(<IN>){
                chomp;
                my $tab="\t"x11;
                my @arr=split (/\t/,$_,3);
                my @gene_list;
                if ( $arr[0] =~ /;/ ) {
                        my @snp_type = split ( /;/, $arr[0]);
                        my @gene_type = split (/;/,$arr[1]);
                        foreach my $index ( 0..$#snp_type ) {
                            my @Gene_List = snp_parse($snp_type[$index],$gene_type[$index]);
                            map { push @gene_list, $_ } @Gene_List;
                        }

                }else{
                        my @Gene_List = snp_parse($arr[0],$arr[1]);
                        map { push @gene_list,$_ } @Gene_List;
                }	
		
		my @uniq = uniq(\@gene_list);
                $gene2ann{$uniq[0]}="\t-\t-\t-\t-\t-\t-\t-" unless $gene2ann{$uniq[0]};
                print SNP_OUT "$_\t$uniq[0]$gene2ann{$uniq[0]}\n";

                if ( @uniq >1 ){
                    foreach my $index ( 1..$#uniq ) {
                        $gene2ann{$uniq[$index]} = "\t-\t-\t-\t-\t-\t-\t-" unless $gene2ann{$uniq[$index]};
                              print SNP_OUT "$tab\t$uniq[$index]$gene2ann{$uniq[$index]}\n";
                    }
                }
          }

	}elsif ( $filename =~ /\.exonic_variant_function\.xls/ ) {
            open IN,$in or die $!;
	    open SNP_OUT,">$outfile" or die $!;
            chomp ( my $head = <IN> );
            print SNP_OUT "$head\t(GeneName)GeneSymbol\tGeneID\tDescription\tBiological Process\tCellular Component\tMolecular Function\tKO_ID(KO_name)\tKEGG_Pathway\n";;
            while(<IN>){
                chomp;
                my @arr = split (/\t/,$_,4);
                my $gene = (split /:/,$arr[2])[0];
                $gene2ann{$gene} = "\t-\t-\t-\t-\t-\t-\t-" unless $gene2ann{$gene};
                #$gene2ann{$gene} = "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-" unless $gene2ann{$gene};
                print SNP_OUT "$_\t$gene$gene2ann{$gene}\n";
        }
   }
	close IN;
	close SNP_OUT;

       ###rename *.ann files###
       if ( $type ) {
       my $dir = dirname $outfile;
       `rename '.ann' '' $dir/*.ann`;
       }
}

sub USAGE{
	my $usage=<<"USAGE";
Program: $0

Contactor: Jiangdezhi (jiangdezhi\@berrygenomics.com)

Usage:
	-indir	<dir>	result directory of RNA_seq_ref (forced)
	-ann	<file>	gene annotation file (forced)
	-SNP		snp analysis flag 
	-AS		altersplice analysis flag 
	
Example:
	perl $0 -indir -ann -SNP -AS

USAGE
	print "$usage\n";
	exit;
}


