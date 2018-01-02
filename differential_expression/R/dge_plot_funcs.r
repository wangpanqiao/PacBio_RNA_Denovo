### load library
require(ggplot2)


#################################################################
### Drawing MA plot for differential expression analysis.
### For each gene: the log2(fold change) between the two samples is plotted (y axis)
### against the gene's log2(average expression) in the two samples(x axis)
### NOTE: here using RPKM stand for expression
#################################################################
plot_MA <- function(log10exp=NULL, log2FC=NULL, FDR=NULL, Significant=NULL,
xlab="log10(RPKM)", ylab="log2(FC)", main="MA plot") {
	# check args
	# check null
	if( is.null(log2FC) ) stop("log2FC is NULL")
	if( is.null(FDR) ) stop("FDR is NULL")
	if( is.null(Significant) ) stop("Significant is NULL")
	# check length
	len <- c(length(log10exp), length(log2FC), length(FDR))
	if( len[1] != len[2] ) 
		stop(paste("length(log10exp) != length(log2FC): ", len[1], " != ", len[2], sep=""))
	if( len[2] != len[3] ) 
		stop(paste("length(log2FC) != length(FDR): ", len[2], " != ", len[3], sep=""))
	if( len[1] == 0 ) stop("length(log2FC) == 0")

	# plot
	p <- qplot(log10exp, log2FC, xlab=xlab, ylab=ylab, main=main, size=I(0.5), colour=Significant)

	# return
	return(p)
}


#################################################################
### Drawing Volcano plot for differential expression analysis.
### For each gene: the log2(fold change) between the two samples is plotted (x axis)
### against the gene's -log10(FDR) in the two samples(y axis)
#################################################################
plot_Volcano <- function(log2FC=NULL, FDR=NULL, Significant=NULL, 
xlab="log2(FC)", ylab="-log10(FDR)", main="Volcano plot") {
	# check args
	# check null
	if( is.null(log2FC) ) stop("log2FC is NULL")
	if( is.null(FDR) ) stop("FDR is NULL")
	if( is.null(Significant) ) stop("Significant is NULL")
	# check length
	len <- c(length(log2FC), length(FDR))
	if( len[1] != len[2] ) 
		stop(paste("length(log2FC) != length(log2FDR): ", len[1], " != ", len[2], sep=""))
	if( len[1] == 0 ) stop("length(log2FC) == 0")

	# plot
	p <- qplot(log2FC, -log10(FDR), xlab=xlab, ylab=ylab, main=main, size=I(0.5), colour=Significant)

	# return
	return(p)
}



#################################################################
### Drawing MA plot and Volcano plot for differential expression analysis.
### NOTE: here using RPKM stand for expression
#################################################################
plot_MA_Volcano <- function(ma_file=NULL, vo_file=NULL, log10exp=NULL, log2FC=NULL, FDR=NULL, Significant=NULL,
exp_lab="log2(RPKM)", FC_lab="log2(FC)", FDR_lab="-1*log10(FDR)", ma_main="MA plot", vo_main="Volcano plot") {
	# check null
	if( is.null(ma_file) ) stop("ma_file is NULL")
	if( is.null(vo_file) ) stop("vo_file is NULL")

	# MA plot
	ma <- plot_MA(log10exp=log10exp, log2FC=log2FC, FDR=FDR, Significant=Significant,
		xlab=exp_lab, ylab=FC_lab, main=ma_main)

	# Volcano plot
	vo <- plot_Volcano(log2FC=log2FC, FDR=FDR, Significant=Significant,
		xlab=FC_lab, ylab=FDR_lab, main=vo_main)

	# output the plot into the file with PNG format
	# for MA plot
	png(filename=ma_file, height = 3000, width = 3400, res = 500, units = "px")
	print(ma)
	dev.off()
	# for Volcano plot
	png(filename=vo_file, height = 3000, width = 3400, res = 500, units = "px")
	print(vo)
	dev.off()

	# return
	return(list(MA=ma, Volcano=vo))
}



#################################################################
### Drawing correlation plot for expression of two samples.
### NOTE: here using RPKM stand for expression
#################################################################
plot_corr <- function(file=NULL, rpkm1=NULL, rpkm2=NULL, xlab="log10(RPKM1)", ylab="log10(RPKM2)", main="correlation plot") {
	# check args
	# check null
	if( is.null(file) ) stop("file is NULL")
	if( is.null(rpkm1) ) stop("rpkm1 is NULL")
	if( is.null(rpkm2) ) stop("rpkm2 is NULL")
	# check length
	len <- c(length(rpkm1), length(rpkm2))
	if( len[1] != len[2] ) 
		stop(paste("length(rpkm1) != length(rpkm2): ", len[1], " != ", len[2], sep=""))
	if( len[1] == 0 ) stop("length(rpkm1) == 0")

	# plot
	p <- qplot(log10(rpkm1), log10(rpkm2), xlab=xlab, ylab=ylab, main=main, size=I(0.5))
	p <- p + geom_abline(slope=1, size=0.5, colour=6, alpha=0.5)
	#p <- p + xlim(-2, 6) + ylim(-2, 6)

	# output the plot into the file with PNG format
	png(filename=file, height = 3000, width = 3000, res = 500, units = "px")
	print(p)
	dev.off()

	# return
	return(p)
}


