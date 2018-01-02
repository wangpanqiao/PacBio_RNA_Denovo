rm(list=ls())
args <- commandArgs(trailingOnly = TRUE )
if ( length(args)<2 ){
                     stop("Please input the information of parameters:
                     args[1]	infile
                     args[2]	outfile
                     args[3]	outdir
		     args[4]	title of graph 
 
                     Usage:
                            Rscript kobas_RNA_single.R infile outfile outdir title
		     Example:
			    Rscript kobas_RNA_single.R /share/work2/staff/jiangdezhi/Program/RNA_seq_ref/BFC2014354_medicago_truncatula/BFC2014354/rh1_b_case_vs_WT_b_ctrl_edgeR/Kegg_Pathway/rh1_b_case_vs_WT_b_ctrl_Up.identify pathway_enrichment.png ./ 'up regulation' " 
        )
}
infile<-args[1]
outfile<-args[2]
outdir<-args[3]
prefix<-args[4]

setwd(outdir)
png(outfile, width=1000, height=800)
kegg<-file(infile, 'r')
rawdata<- readLines(kegg)
data<-rawdata[-which(grepl("^(#|-|\\/)|^$",rawdata),arr.ind<-TRUE)]
#data<-rawdata[-which(grepl("^(#|-)|^$",rawdata),arr.ind<-TRUE)]
if ( length(data) >20 ){
      num<-1:20
      data_list<-sapply(data[num],function(X){strsplit(X,"\t")})
}else{
      data_list<-sapply(data,function(X){strsplit(X,"\t")})
}

term<-vector()
gene_num<-vector()
for (i in 1:length(data_list)){
                         term[i]<-data_list[[i]][1]
                         gene_num[i]<-as.numeric(data_list[[i]][4])
}

term_gene_num=as.matrix(gene_num)
par(mar=c(15,10,3,5),mgp=c(3,1,0))
position <- barplot(term_gene_num[,1], 
                    ylim=c(0, 1.2*max(term_gene_num)),
                    ylab="Gene Number", 
                    cex.lab=1.5,
                    main=prefix,
                    col="#6B8E23"
)

text(position, -0.01,
     labels =term,adj=c(1.0,1.0),
     srt = 30,
     cex=1.2,
     xpd = NA,
)
text(position,
     term_gene_num[,1],
     adj=0.5,
     labels=term_gene_num[,1],
     pos=3
)
box()
dev.off()
