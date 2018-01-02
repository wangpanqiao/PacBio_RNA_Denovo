#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use 5.010;
use Data::Dumper;
use DBI;
use File::Spec;
use lib "/home/wudi/lib";
require "ABSyS.pl";

my ($config_file,
	$out_path,
	$prefix);
GetOptions(
	"cfg:s"   => \$config_file,
	"o:s"	  => \$out_path,
	"prefix:s"=> \$prefix,
	"h|?" => sub{
say "Usage:
	perl $0 -cfg < config file > -o < /output/path > -prefix < prefix >
";
	exit(0)
	}
);


my $shell = "$out_path/$prefix.sh";

if( defined $out_path && -d $out_path){
	$out_path = ABSOLUTE_DIR($out_path);
}elsif(defined $out_path){
	MKDIR($out_path);
	$out_path = ABSOLUTE_DIR($out_path);
}
MKDIR("$out_path/ALL_TO_ALL_BLAST");
my ($sample_list,$samples) ;
open SL, "$config_file" || die "Can't find the sample list\n";

####### SOFTWARE WHERE ############

my $FastaU = "/home/wudi/script/Fasta.Utili.pl ";
my $AdjustFasta = " /home/wudi/workdir/KG/TreeFam/bin/orthomclAdjustFasta ";
my $dbformat = "/share/software/software/blast-2.2.20/bin/formatdb ";
my $blastp = "/share/work2/staff/zhouy/Software/Blast/bin/blastp "; 
my $INSTALL_SCHEMA = "/home/zhouy/work2/Software/OrthoMCL/bin/orthomclInstallSchema";
my $BLAST_PARSER = "/home/zhouy/work2/Software/OrthoMCL/bin/orthomclBlastParser";
my $LOAD_BLAST = "/home/zhouy/work2/Software/OrthoMCL/bin/orthomclLoadBlast";
my $PAIRS = "/home/zhouy/work2/Software/OrthoMCL/bin/orthomclPairs";
my $output_format = "\'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \' " ;
my $DUMP_PAIRS = "/home/zhouy/work2/Software/OrthoMCL/bin/orthomclDumpPairsFiles";
my $MCL = "/home/zhouy/work2/Software/MCL/bin/mcl";
my $MCL_GROUP = "/home/zhouy/work2/Software/OrthoMCL/bin/orthomclMclToGroups" ;
my $frm = "/home/wudi/workdir/KG/TreeFam/bin/format_output.pl ";
my $re4R = "/home/wudi/workdir/KG/TreeFam/bin/reformat_4R.pl ";
my $ploter = "/home/wudi/workdir/KG/TreeFam/bin/Protein_Fam.plot.R ";

####### SOFTWARE END ############

####### LOAD sample list #######

my $i = 0;
while(<SL>){
	chomp;
	my @array = split "\t",$_;
	$samples->[$i]->{"tax"} = $array[0];
	$samples->[$i]->{"fasta_file"} = $array[1];
	$samples->[$i]->{"accession_region"} = $array[2];
	$i++ ;
}
close SL;

###### END #####


## MAKE MYSQL CONFIG ## 
if ( create_mysql_config($prefix,"$out_path/mysql.config") ||
    SETUP_DATABASES($prefix) ){
	warn "Some mistake happen in database setup step\n";
}
##### ADjust Fasta all fasta file #####

$i = 0;
foreach ( @{$samples} ){
	my $sample = $_;
	system("$AdjustFasta $sample->{'tax'} $sample->{'fasta_file'} $sample->{'accession_region'}  $out_path/ALL_TO_ALL_BLAST/"); 
	$samples->[$i]->{'fasta_file'} = "$out_path/ALL_TO_ALL_BLAST/$sample->{'tax'}.fasta";
	$i++;
}
$i = 0;

open SHELL, ">$shell";
system("cat ".(join " ",map {$_->{'fasta_file'}} @$samples )." > "."$out_path/$prefix.all.fasta");
system($dbformat." -i "." $out_path/$prefix.all.fasta");

##### make blast sh #####

foreach ( @{$samples} ){
	my $sample = $_;
	open BLAST,">$out_path/ALL_TO_ALL_BLAST/blast.$sample->{'tax'}.sh"; 
	say BLAST $blastp." -evalue 1e-10 -query $sample->{'fasta_file'} ".
	" -out $out_path/ALL_TO_ALL_BLAST/$sample->{'tax'}.out -db".
	" $out_path/$prefix.all.fasta -max_target_seqs 10 -num_threads 14 -outfmt $output_format"; ### FORMAT NEED TO BE FIXED ####
	close BLAST;
	$samples->[$i]->{'blast_result'} = "$out_path/ALL_TO_ALL_BLAST/$sample->{'tax'}.out";
	$i ++;
}

say SHELL "cat ", (join " ",map {$_->{'blast_result'}} @$samples ) ,">" , " $out_path/$prefix.all.blast_result";
say SHELL "mkdir $out_path/ALL_TO_ALL_BLAST/Fasta && mv $out_path/ALL_TO_ALL_BLAST/*.fasta $out_path/ALL_TO_ALL_BLAST/Fasta ";

##### WHAT SHELL DO ######

say SHELL "$INSTALL_SCHEMA $out_path/mysql.config  $out_path/orthomcl.log ";
say SHELL "$BLAST_PARSER $out_path/$prefix.all.blast_result $out_path/ALL_TO_ALL_BLAST/Fasta >> $out_path/similarSequences.txt ";
say SHELL "$LOAD_BLAST $out_path/mysql.config $out_path/similarSequences.txt ";
say SHELL "$PAIRS $out_path/mysql.config $out_path/orthomcl.log cleanup=no ";
say SHELL "$DUMP_PAIRS $out_path/mysql.config "; ## need fixed ##
say SHELL "$MCL mclInput --abc -o $out_path/result.out -I 2.0";
say SHELL "cat $out_path/result.out | $MCL_GROUP GROUP 1 > $out_path/$prefix.temp.group.xls";
### DROP DATA SCRIPT THERE ###

say SHELL "$frm -i $out_path/$prefix.temp.group.xls -o $out_path/$prefix.group.xls -t orthomcl" ;
say SHELL "$re4R ",(join ",",map {$_->{'tax'}} @$samples )," $out_path/$prefix.group.xls"," $out_path/$prefix.4r";
say SHELL "$ploter -i  $out_path/$prefix.4r -n ",(join ",",map {$_->{'tax'}} @$samples ), " -o $out_path/$prefix.venn.diagram.png";
Make();

say "Work Done\n";
exit(0);
#########################################################3

sub Make{
	open MAKE,">$out_path/orthomcl.make";
	say MAKE "orthomcl.end :  $prefix.venn.diagram.png ";
	say MAKE "	touch orthomcl.end \n";
	foreach ( @{$samples} ){
		my $sample = $_;
		say MAKE "ALL_TO_ALL_BLAST/$sample->{'tax'}.out : ";
		say MAKE "	cd ALL_TO_ALL_BLAST/ &&	qsub  -l  vf=3G  -q  all.q -cwd -sync y blast.$sample->{'tax'}.sh	\n";
		$i ++;
	}
	say MAKE "$prefix.venn.diagram.png : ",(join " ",map{ File::Spec->abs2rel($_->{'blast_result'},$out_path) } @$samples );
	say MAKE "	sh ".File::Spec->abs2rel($shell,$out_path);
	say MAKE "clean:";
	say MAKE "	find ./ -name '*.sh.[eo][0-9]*' | xargs rm ";

}

sub SETUP_DATABASES{
	my $sample_name = shift;
	my $dbh = DBI->connect("DBI:mysql:host=10.255.255.246;port=3306;mysql_local_infile=1",
							'pacbio',
							"pacbio123") || die "Can't Connect host\n";

my $SQL1=<<SQL1;
drop database if exists \`$sample_name\`;
SQL1
my $SQL2=<<SQL2;
create database $sample_name; 
SQL2
	my $sth = $dbh->prepare($SQL1);
    my $sta = $sth->execute() || die "Can't Drop the previous $sample_name database";
	$sth = $dbh->prepare($SQL2);
	$sta = $sth->execute() || die "Can't Create the $sample_name database";

	$dbh->disconnect();
	return 0;
}

sub create_mysql_config{

	my $sample_name = shift;
	my $output_file = shift;
	my $TEMPLATE=<<TEMP;
# this config assumes a mysql database named 'orthomcl'.adjust according
# to your situation.
dbVendor=mysql 
dbConnectString=dbi:mysql:$prefix:10.255.255.246:3306;mysql_local_infile=1
dbLogin=pacbio
dbPassword=pacbio123
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE
TEMP
	open FH,">$output_file" || die "Can't create the mysql config file\n";
	print FH $TEMPLATE;
	close FH;
	return 0;
}

__END__
