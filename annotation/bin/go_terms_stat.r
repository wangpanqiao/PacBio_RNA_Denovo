# usage function
usage <- function() {
        print("-------------------------------------------------------------------------------")
        print("Usage: Rscript go_terms_stat.r gene2go.map out_prefix")
        print("1) gene2go.map: the genes go annotation list")
        print("2) out_prefix: the png & pdf files prefix")
        print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) < 2 ) {
        print(args)
        usage()
        stop("the length of args < 2")
}

infiles <- args[1]
outfile <- args[2]

.libPaths("/home/jiangdezhi/bin/R_installed_package/")
library(topGO)
geneID2GO <- readMappings(infiles)
str(head(geneID2GO))

geneNames <- names(geneID2GO)
myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

TopGO_ontology<-function(GeneList,ontology="BP",Outfile=outfile){
	GOdata<-new("topGOdata",description=paste('allgenes',ontology,sep='_'),ontology=ontology, allGenes=GeneList,annot = annFUN.gene2GO,gene2GO = geneID2GO)
	GOstat<-termStat(GOdata)
	write.table(GOstat,file=paste(Outfile,ontology,'TermStat.xls',sep='_'),col.names=TRUE,append=TRUE,quote=F,sep='\t')

	return(GOstat)
}

BPdata<-TopGO_ontology(geneList,"BP",outfile)
MFdata<-TopGO_ontology(geneList,"MF",outfile)
CCdata<-TopGO_ontology(geneList,"CC",outfile)
