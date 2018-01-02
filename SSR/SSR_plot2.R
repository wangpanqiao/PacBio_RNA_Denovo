#===================================================================
args<-commandArgs(TRUE)
input=args[1]
dir=args[2]
library(ggplot2)
library(reshape2)
library(gplots)
library(scatterplot3d)
library(RColorBrewer)
t<-read.delim(input,sep='\t', header=F,stringsAsFactors=FALSE)
if(dim(t)[2] > 12){
    t <- t[,1:12]
}
colnames(t) = t[1,]
t = as.data.frame(t[2:nrow(t),])
type_columns = 2:ncol(t)
colnames(t)[type_columns] = gsub("repeat\\((.*)\\)", "\\1", colnames(t)[type_columns])
t[type_columns] = lapply(t[type_columns], as.numeric)
is.zero = which(apply(t[type_columns], 2, sum) == 0) + 1
if(length(is.zero)==0){
}else{
    t = t[,-is.zero]
}
df = melt(t, id.vars="number")
df$number = factor(df$number, levels= t$number)
names(df)[1] = "repeat_typerweg"
name=as.character(t[[1]])
length =length(names(t)[-1])
if(length<=8){
    my_color1 <- (brewer.pal(length,"Set2"))
}else{
    my_color1 <- (brewer.pal(length,"Set3"))
}
png(file= paste(dir, 'Transcripts.fa.SSR_motifs_distribution.png',sep='/'), type="cairo-png")
scatterplot3d(df,type="h",lwd=12,pch="",x.ticklabs=name,xlab=("SSR motif unit"),ylab=("repeat_type"),zlab=("Repeat counts"),box = FALSE,color=(rep(my_color1, each=length(name))),main="Distribution of SSR Motifs",y.margin.add=0.5,angle=23)
lengname<-names(t)[-1]
legend("topright",pch=c(rep(16,length)),legend=lengname,col=my_color1)
pdf(file= paste(dir, 'Transcripts.fa.SSR_motifs_distribution.pdf',sep='/'))
scatterplot3d(df,type="h",lwd=12,pch="",x.ticklabs=name,xlab=("SSR motif unit"),ylab=("repeat_type"),zlab=("Repeat counts"),box = FALSE,color=(rep(my_color1, each=length(name))),main="Distribution of SSR Motifs",y.margin.add=0.5,angle=23)
legend("topright",pch=c(rep(16,length)),legend=lengname,col=my_color1)
#===================================================================
