#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
#####################################################################
#	This program is used for plotting SSR class.
#			Bar plot.
#
#####################################################################
die "Usaage: perl $0 <ssr_classification.txt>\nThe figures and plotting R script are in current dir\n" unless @ARGV==1;

my $stat=shift;

#my $f1=0;

open STAT,$stat;
my $first=<STAT>;
my @array=split /\t/,$first;
my $size=$array[1];
#open OUT,">$name.SSR_density.txt";
#while(<STAT>){
#	chomp;
#	if(/^Total size of examined sequences \(bp\)\:\s+(\d+)/){
#		$size=$1;
#	}
#	if(/^Unit/){
#		$f1=1;
#	}elsif(/^\s+/ && $f1==1){
#		last;
#	}
#
#	if($f1==1){
#		print OUT $_;
#	}
#}
#close OUT; 
close STAT;

my $R= <<R;
#============================================================================================================
dat<-read.table("$name.SSR_density.txt",header=TRUE,row.names=1,sep="\\t")
dat2<-round(dat*1000000/$size);
dat2<-as.matrix(dat2)
colnames(dat2)<-NULL;
png("$name.SSR_density.png",type="cairo-png")
barplot(dat2,beside=T,ylab="SSR count/bases (MB)",col="#983063",space=0.75,ylim=c(0,max(dat2)*1.2),main="SSR density")
axis(side=1,at=(0:nrow(dat2)*1.75+0.5),labels=rep("",nrow(dat2)+1))
text(1:nrow(dat2)*1.75-0.5,-round(max(dat2)*1.2)/25,labels=rownames(dat2),xpd=T)
mtext(side=1,line=2,"Period number")
dev.off()


pdf("$name.SSR_density.pdf")
barplot(dat2,beside=T,ylab="SSR count/bases (MB)",col="#983063",space=0.75,ylim=c(0,max(dat2)*1.2),main="SSR density")
axis(side=1,at=(0:nrow(dat2)*1.75+0.5),labels=rep("",nrow(dat2)+1))
text(1:nrow(dat2)*1.75-0.5,-round(max(dat2)*1.2)/25,labels=rownames(dat2),xpd=T)
mtext(side=1,line=2,"Period number")
dev.off()


#============================================================================================================
R

open R,"|/PUBLIC/software/public/System/R-2.15.3/bin/R --slave --vanilla" or die $!;
print R $R;
close R;

open R,">ssr_plot.R";
print R $R;
close R;
