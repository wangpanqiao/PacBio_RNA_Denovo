#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
#####################################################################
#	This program is used for plotting SSR class.
#			Bar plot.
#
#####################################################################
die "Usaage: perl $0 <SSR.statistics>\nThe figures and plotting R script are in current dir\n" unless @ARGV==1;

my $stat=shift;
my $name=basename($stat,".statistics");
my $size;

my $f1=0;

open STAT,$stat;
open OUT,">SSR_density.txt";
while(<STAT>){
#	chomp;
	if(/^Total size of examined sequences \(bp\)\:\s+(\d+)/){
		$size=$1;
	}
	if(/^Unit/){
		$f1=1;
	}elsif(/^\s+/ && $f1==1){
		last;
	}

	if($f1==1){
		print OUT $_;
	}
}
close OUT; close STAT;

my $R= <<R;
#============================================================================================================
dat<-read.table("$name.SSR_density.txt",header=TRUE,row.names=1,sep="\\t")
dat2<-round(dat*1000000/$size);
dat2<-as.matrix(dat2)
colnames(dat2)<-NULL;
png("SSR_density.png",type="cairo-png")
par(mar=c(4, 6, 4, 2))
barplot(dat2,beside=T,ylab="SSR count/bases (MB)",cex.lab=1.8,col="#983063",space=0.75,ylim=c(0,max(dat2)*1.2),main="SSR density", font=1,cex.axis=1.5,cex.main=2)
axis(side=1,at=(0:nrow(dat2)*1.75+0.37),labels=rep("",nrow(dat2)+1))
text(1:nrow(dat2)*1.75-0.5,-round(max(dat2)*1.2)/25,labels=rownames(dat2),xpd=T,cex=1.5)
mtext(side=1,line=2.5,"Period number",cex=1.8)
dev.off()


pdf("SSR_density.pdf")
par(mar=c(4, 6, 4, 2))
barplot(dat2,beside=T,ylab="SSR count/bases (MB)",cex.lab=1.8,col="#983063",space=0.75,ylim=c(0,max(dat2)*1.2),main="SSR density", font=1,cex.axis=1.5,cex.main=2)
axis(side=1,at=(0:nrow(dat2)*1.75+0.37),labels=rep("",nrow(dat2)+1))
text(1:nrow(dat2)*1.75-0.5,-round(max(dat2)*1.2)/25,labels=rownames(dat2),xpd=T,cex=1.5)
mtext(side=1,line=2.5,"Period number",cex=1.8)
dev.off()






#============================================================================================================
R

open R,"|/share/public/software/R-2.15.2/bin/R --slave --vanilla" or die $!;
print R $R;
close R;

open R,">ssr_plot.R";
print R $R;
close R;
