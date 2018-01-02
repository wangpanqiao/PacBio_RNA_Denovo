#!/usr/bin/perl
#
use File::Basename;
die "Usage: perl $0 <ssr_motif_frequency.txt> \n" unless @ARGV==1;

my $frequency=shift;

my $dir=dirname($frequency);
my $name=basename($frequency,".ssr_motif_frequency.txt");
open(SSR,"$frequency");

my @start,@end,@start_pos,@end_pos;
my $head=<SSR>;
chomp $head;
my @head=split /\t/,$head;
my $arrar_num=4;
for ($i=1;$i+$arrar_num-1<$#head;$i+=$arrar_num){
#	print $i+1,"\t",$head[$i],"\t",$i+$arrar_num,"\t",$head[$i+$arrar_num-1],"\n";
	push @start,$head[$i];
	push @start_pos,$i;
	push @end,$head[$i+$arrar_num-1];
	push @end_pos,$i+$arrar_num-1;
}
$end[$#end]=$head[$#head-1];
$end_pos[$#end_pos]=$#head-1;

my @repeat_name;
my @lengend_name;
for my $i(0..$#start){
	my $tmp="repeat(".$start[$i]."-".$end[$i].")";
	my $leng=$start[$i]."-".$end[$i];
	push @repeat_name,$tmp;
	push @lengend_name,$leng;
}
my $lengend_name=join ",",@lengend_name;
$lengend_name="\"".$lengend_name."\"";

my %hash;
while(<SSR>){
	chomp;
	my @tmp=split /\t/, $_;
	my $num=length($tmp[0]);
	for my $i(0..$#repeat_name){
		for my $j($start_pos[$i]..$end_pos[$i]){
			if(defined $tmp[$j]){$hash{$num}{$repeat_name[$i]}+=$tmp[$j];}
			else{$hash{$num}{$repeat_name[$i]}+=0;}
		}
	}
}
close SSR;

my @arrra=("Mono-","Di-","Tri-","Tetra-","Penta-","Hexa-","Hepta-","Octo-","Nona-","Deca-");
open(OUT,">$dir/SSR_summary_plot.txt");
print OUT "number\t",join("\t",@repeat_name),"\n";
foreach my $num(sort keys %hash){
	my @tmp;
	for my $i(0..$#repeat_name){
		push @tmp,$hash{$num}{$repeat_name[$i]};
	}
	print OUT $arrra[$num-1],"\t",join("\t",@tmp),"\n";
}

close OUT;
my $R= <<R;
#===================================================================
library(ggplot2)
library(reshape2)
setwd("$dir")
t<-read.delim("SSR_summary_plot.txt",sep='\\t')
t<-as.data.frame(t)
colnames(t)<-c("repeat_type",paste(1:(length(t)-1),colnames(t)[2:length(t)],sep=":"))
df<-melt(t)
names_rep_leg<-unlist(strsplit($lengend_name,","))
name=as.character(t[[1]])
p<-ggplot(df,aes(x=repeat_type,y=value))+geom_bar(aes(fill=variable),width=0.7)+theme(axis.text.x=element_text(colour="black",size=8),axis.text.y=element_text(size=8),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),plot.title=element_text(size=12))+ggtitle("Distribution of SSR Motifs\n")+xlab("SSR motif unit")+ylab("Repeat counts\n")+scale_fill_brewer("repeat_type",label=names_rep_leg,palette="Set2")+xlim(name)
#p<-ggplot(df,aes(x=nucl_num,fill=repeat_type))
#p1<-p+geom_bar(width=0.7)+theme(axis.text.x=element_text(colour="black",size=8),axis.text.y=element_text(size=8),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),plot.title=element_text(size=12))+xlim(names_nucl)+ggtitle("Distribution of SSR Motifs\\n")+xlab("SSR motif unit")+ylab("Repeat counts\\n")+scale_fill_brewer(label=names_rep_leg,palette="Set2")
ggsave("$name.SSR_motifs_distribution.png",type="cairo-png",width=6,height=4,plot=p)
ggsave("$name.SSR_motifs_distribution.pdf",width=6,height=4,plot=p)

#===================================================================
R

open R,"|/PUBLIC/software/public/System/R-2.15.3/bin/R --slave --vanilla"  or die $!;
print R $R;
close R;

open R,">SSR_plot2.R";
print R $R;
close R;
