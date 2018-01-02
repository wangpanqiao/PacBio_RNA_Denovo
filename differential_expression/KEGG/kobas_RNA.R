rm(list=ls())
args<-commandArgs(trailingOnly = T)
if(length(args)<3) stop ("Uasge:
	Rscript kobas_new.R Up_pathway_enchriment.xls Down_pathway_enchriment.xls outfile"
)

up<-args[1]
down<-args[2]
outfile<-args[3]

data_up<-read.delim(up,head=T,sep="\t",stringsAsFactors=F)
data_down<-read.delim(down,head=T,sep="\t",stringsAsFactors=F)
			
if(length(data_up[,4])>20){ 
	up<-data_up[,4][1:20] 
}else{ 
	up<-data_up[,4]
}

if(length(data_down[,4])>20){
	down<-data_down[,4][1:20]
}else{
	down<-data_down[,4]
}

enriched_pathway <- c(up,down)			 
max<-max(as.numeric(enriched_pathway))

png(outfile,height=1000,width=1600)
par(mfrow=c(1,2), mar=c(22,10,2,0),las=1, cex=2)

layout(matrix(c(1,2), nrow=1), width=c(length(up), length(down)))
up_position<- barplot(as.numeric(up), ylim=c(0, 1.1*max),col="#6B8E23", main='Up regulation')
mtext("Gene Count", at=1.1*max/2, side=2, line=3, las=0, cex=1.3)
text(up_position, -0.1, labels=data_up[,1][1:length(up)], srt=45, adj=1, xpd=NA, cex=1.2)
text(up_position, as.numeric(up), labels=up, pos=3)
box()

par(mar=c(22,0,2,10))
down_position <- barplot(as.numeric(down), ylim=c(0, 1.1*max), col="#A0522D",yaxt='n', main='Down regulation')
text(down_position, -0.1, labels=data_down[,1][1:length(down)], srt=45, adj=1, xpd=NA, cex=1.2)
text(down_position, as.numeric(down), labels=down, pos=3)
box()
dev.off()
