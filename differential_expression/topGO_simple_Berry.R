rm(list=ls())
getProgramName<-function(){
        args <- commandArgs(trailingOnly = FALSE)
        sub("--file=", "", args[grep("--file=", args)])
}
program <- getProgramName()
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<2) {
	stop(sprintf("Please input the information of arguments:
	args[1]	<file>	infile
	args[2]	<file>	gene2GO
	args[3]	<str>	outfile_prefix
	args[4]	<dir>	outdir
        
Example:
	Rscript %s /share/work2/staff/jiangdezhi/Program/RNA_seq_ref/BFC2014143_rat/BFC2014143/Specific_Analysis/M_Y_G_Z_D_Analysis/diff_gene.xls /share/work2/staff/jiangdezhi/Program/RNA_seq_ref/BFC2014143_rat/BFC2014143/Specific_Analysis/M_Y_G_Z_D_Analysis/gene2go.xls topGO ./
Usage:
	Rscript %s infile gene2GO outfile_prefix outdir",program,program)
	)
}

ptm <- proc.time()

infile<- args[1]
soybean_gene2GO_infile<-args[2]
outfile<- args[3]
outdir<- args[4]

###load topGO package
#.libPaths("/home/xuxiong/bin/R_installed_package/")
#library(topGO)
#library(AnnotationFuncs)
.libPaths("/home/jiangdezhi/bin/R_installed_package/")
library(topGO)
library(ggplot2)


###load gene2go and input file
soybean_gene2GO<-read.delim(soybean_gene2GO_infile,header=F,sep="\t",comment.char = '#',blank.lines.skip=TRUE, encoding = "UTF-8",stringsAsFactors=FALSE)
rawData<-read.table(infile,header=T,sep="\t",comment.char = '#',blank.lines.skip=TRUE, encoding = "UTF-8",stringsAsFactors=FALSE)

###create geneList for GO analysis
genes<-rawData[,1]
soybean_gene2GO_2_10<-unique(soybean_gene2GO[,c(1,2)])
L_geneID2GO<-strsplit(soybean_gene2GO_2_10[,2],", ")
names(L_geneID2GO)<-soybean_gene2GO_2_10[,1]
geneList<-rep(1,length(L_geneID2GO))
geneList[names(L_geneID2GO) %in% rawData[,1]]=0.05
names(geneList)<-names(L_geneID2GO)

###sub function
topDiffGenes<-function(allScore){return(allScore<=0.05)}#select differential genes

TopGO_ontology<-function(GeneList,ontology="BP",Outfile=outfile){#create object for GO analysis
	GOdata<-new("topGOdata",description=paste('diffgenes',ontology,sep='_'),ontology=ontology,geneSel=topDiffGenes,allGenes=GeneList,annot = annFUN.gene2GO,gene2GO = L_geneID2GO)
	resultFisher<-runTest(GOdata,algorithm="classic",statistic="fisher")
	resultKS<-runTest(GOdata,algorithm="classic",statistic="ks")
	resultKS.elim<-runTest(GOdata,algorithm="elim",statistic="ks")
	allRes<-GenTable(GOdata,classicFisher=resultFisher,classicKS=resultKS,elimKS=resultKS.elim,orderBy="elimKS",ranksOf="classicFisher",topNodes=attributes(resultKS.elim)$geneData[4],numChar=80)
	allRes$Term<-gsub("\\s+",'_',allRes$Term,perl=TRUE)
	ann.genes_list<-genesInTerm(GOdata, allRes$GO.ID)
	sig.genes<-sigGenes(GOdata)
	allRes$ann.genes <- sapply(allRes$GO.ID,function(X){
			sig<-ann.genes_list[[X]] %in% sig.genes
                        paste(ann.genes_list[[X]][sig], collapse = ';')
			})
	print(allRes$ann.genes)
	select.allRes<-allRes[,c("GO.ID","Term","Annotated","Significant","elimKS","ann.genes")]
        colnames(select.allRes)<-c("GO.ID","Term","Annotated","Significant","PValue","ann.genes")
	cat(c(ontology,"\t"),file=paste(Outfile,ontology,'GenTable.xls',sep='_'))
	write.table(select.allRes,file=paste(Outfile,ontology,'GenTable.xls',sep='_'),col.names=TRUE,append=TRUE,quote=F,sep='\t')

	#png(file =paste(Outfile,ontology,'pValue.classic_VS_pValue.elim.png',sep='_'),res=150,height=600,width=600)
	pValue.classic<-score(resultKS)
	pValue.elim<-score(resultKS.elim)[names(pValue.classic)]
	gstat<-termStat(GOdata,names(pValue.classic))
	gSize<-gstat$Annotated/max(gstat$Annotated)*4
	gCol<-rainbow(length(gstat$Significant),start=0,end=1/6)
	#plot(pValue.classic,pValue.elim,xlab="p-valueclassic",ylab="p-valueelim",pch=19,cex=gSize,col=gCol)
	#dev.off()

	png(file =paste(Outfile,ontology,'SigOfNodes.png',sep='_'),res=150,height=900,width=900)
	showSigOfNodes(GOdata,score(resultKS.elim),firstSigNodes=5,useInfo= 'all')
	pic=dev.off()
	pdf(file =paste(Outfile,ontology,'SigOfNodes.pdf',sep='_'),height=9,width=9,pointsize=40)
	par(ps=40)
	showSigOfNodes(GOdata,score(resultKS.elim),firstSigNodes=5,useInfo= 'all')
	pic=dev.off()
	return(allRes)
}

BP_CC_MF_barplot<-function(GeneCount,TERM){
	colors <- c('#4682B4','#87CEEB','#6B8E23','#A0522D','#FF8C00','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	par(mar = c(14,12,4.1,3.1))
	mp <- barplot(
		as.numeric(GeneCount),
	#	width=1.2,
	#	legend= colnames(GeneCount),
	#	names.arg= colnames(GeneCount),
	#	space=c(0.1,rep(c(rep(0.1,nrow(GeneCount)-1),1),nrow(GeneCount)-1),rep(0.1,nrow(GeneCount)-1)),
		beside=TRUE,#是否柱状交错
	#	xlab="Percentage_of_mapped_reads(%)",
		ylab="Gene Count",
		col=rep(colors[1:3],each=nrow(GeneCount),len=nrow(GeneCount)*ncol(GeneCount)),
		axisnames=TRUE,
		axes=TRUE,
		plot=TRUE,
		xpd=FALSE,#是否超出
		ylim=c(0,max(as.numeric(GeneCount))*1.2),
	#	cex.lab=1,#伸展坐标比例
		pch=15,
	)
	legend("topleft",#"bottomright",
		legend=colnames(GeneCount),
		cex=1,
		fill=colors[1:3],
		col=colors[1:3],
		inset = 0.01
	)
	abline(v = mp[seq(from=nrow(GeneCount),to=length(mp),by=nrow(GeneCount))]+0.5, col = "black", lwd = 1,lty=4)
	title(main=list("Gene Ontology Analysis"));
	text(mp,as.vector(as.numeric(GeneCount)),labels = as.vector(as.numeric(GeneCount)),adj=c(0.5,-0.5), cex=0.8,xpd = TRUE)
	#print(mp)
	text(mp, par("usr")[3], labels = as.vector(TERM), srt = 35,cex=0.7, adj = c(1,1), xpd = TRUE)
	box();
	pic=dev.off()
}

###result
if (!file.exists(outdir)) {
        dir.create(path=outdir)
}
#GO enrichment analysit
setwd(outdir)
BPdata<-TopGO_ontology(geneList,"BP",outfile)
MFdata<-TopGO_ontology(geneList,"MF",outfile)
CCdata<-TopGO_ontology(geneList,"CC",outfile)

#create GO enrichment barplot graph
L_top20<-list(Biological_Process=BPdata[1:20,],Molecular_Function=MFdata[1:20,],Cellular_Component=CCdata[1:20,])
GeneCount<-sapply(L_top20,function(X){X$Significant})
TERM<-sapply(L_top20,function(X){X$Term})
png(paste(outfile,"GO_Term.png",sep='_'),pointsize=18,width=1600,height=900)
BP_CC_MF_barplot(GeneCount,TERM)
pdf(paste(outfile,"GO_Term.pdf",sep='_'),width=16,height=9)
BP_CC_MF_barplot(GeneCount,TERM)

print(proc.time() - ptm)
