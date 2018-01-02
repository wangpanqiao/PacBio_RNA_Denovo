#========================================================================================================
#args[1]=GO_classification_count.txt
##args[2]=output dir
args<-commandArgs(TRUE)
setwd(args[2])
library(RColorBrewer)
dat<-read.table(args[1],sep="\t",skip=1)
dat = dat[order(dat[,3],dat[,2]),]
#labels<-as.character(count$name)

## * ############################### ggplot2 ###############################
.libPaths("/share/software/software/R-3.1.0_installdir/lib64/R/library");
library(ggplot2)
library(reshape2)
df = dat[,-1]
names(df) <- c("label","GO","number")
df = melt(df,id.vars=c("label","GO"))
df <- transform(df,label=paste(label,"  ",sep=""))
ymax = max(df$value)*1.2
p<-ggplot(df, aes(x=label,fill=GO,y=value)) + geom_bar(stat='identity',position="dodge") +
    theme_classic() + theme(axis.text.x=element_text(angle=60,hjust=1,size=7.5)) + facet_grid(.~GO,shrink=F,scales='free',space='free') +
    xlab("") + ylab("Number of genes")+
	ggtitle("Gene Function Classification (GO)
")+coord_cartesian(ylim=c(0,ymax*1.1)) + scale_y_continuous(breaks=pretty(1:ymax))+theme(legend.position="none")


ggsave("GO_classification_bar.png",type="cairo-png",width=12,height=6,plot=p)
ggsave("GO_classification_bar.pdf",width=12,height=6,plot=p)

#========================================================================================================
