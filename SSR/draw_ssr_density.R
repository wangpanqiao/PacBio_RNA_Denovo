###############################################################################
# Author: Daqing Ai
# Email:  aidaqing@berrygenomics.com
# Last modified : 2014-07-17 16:51
# Filename   : ssrDensityHistogram.R
# Description    : 
###############################################################################

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
	stop("
Example:
	Rscript draw_ssr_density.R ssr_density.txt
	")
}

rnaWide <- read.delim(args[1], header=FALSE)
rnaWide <- t(rnaWide)

chartname <- sub('.*?([^/]+)\\.\\w+$', '\\1.png', args[1])
#cat("Charting the SSR density distribution on ", getwd(), "/", chartname, "\n", sep="")

png(filename=chartname, width=500, height=500)
par(mar=c(5, 5, 3, 0))
snpbar <- barplot(as.numeric(rnaWide[2,]), names.arg=rnaWide[1,], las=1, space=0.5, col='steelblue', ylim=c(0, 1.1*max(as.numeric(rnaWide[2,]))), xlab="Period Number", ylab="SSR count/bases (MB)", main="SSR Density", cex=1.2, cex.axis=1.2, cex.lab=1.5, cex.main=1.5)
text(snpbar, as.numeric(rnaWide[2,]) + 1, labels=rnaWide[2,])
box()
dev.off()
