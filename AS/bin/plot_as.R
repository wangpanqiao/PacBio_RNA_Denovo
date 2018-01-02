#!/share/software/software/R-3.2.0_INSTALL/bin/Rscript
remove(list=ls())
library(optparse)

option_list <- list(
        make_option(c("-i", "--in"), action="store", type="character", help="The list of read alignments"),
        make_option(c("-o", "--outdir"), action="store", type="character", help="The outdir of results")
)

parser <- OptionParser(option_list = option_list, description="Author: chaijingchao561@berrygenomics.com", epilogue="Example: Rscript  %prog -i 6.1step.clean.sortS.uniqS.sortQ.uniqQ.txt -o ./")
opt <- parse_args(parser, positional_arguments = F)

data <- read.table(opt$i, head = T, sep = "\t")
n <- dim(data)[1]
for(i in 1:n){
	png(paste(opt$o, "/", data[i,1], "_", data[i,2], "_as.png", sep = ""), width = 1200, height = 800)
	#plot(c(0,max(data[i,4:11])), c(0, 40), type = "n", xlab = "", ylab = "", axes = F)
	if(data[i,9] < data[i,10]){
			gap = data[i,10] - data[i,9]
			plot(c(0,max(data[i,4:11], data[i,6] + gap, data[i,7] + gap)), c(0, 40), type = "n", xlab = "", ylab = "", axes = F)
			rect(data[i,4], 28, data[i,5], 32, col = "blue")
			rect(data[i,6] + gap, 28, data[i,7] + gap, 32, col = "blue")
			lines(c(data[i,5], data[i,6] + gap), c(30,30))
			lines(c(data[i,4], data[i,8]), c(28, 12), lty = 2)
			lines(c(data[i,5], data[i,9]), c(28, 12), lty = 2)
			lines(c(data[i,6] + gap, data[i,10]), c(28, 12), lty = 2)
			lines(c(data[i,7] + gap, data[i,11]), c(28, 12), lty = 2)
			text(c(data[i,4], data[i,5], data[i,6] + gap, data[i,7] + gap), c(33, 33, 33, 33), c(data[i,4], data[i,5], data[i,6], data[i,7]))
			text(c(data[i,8], data[i,9], data[i,10], data[i,11]), c(7, 7, 7, 7), c(data[i,8], data[i,9], data[i,10], data[i,11]))

	}else{
			gap = data[i,8] - data[i,11]
			plot(c(0,max(data[i,4:11], data[i,4] + gap, data[i,5] + gap)), c(0, 40), type = "n", xlab = "", ylab = "", axes = F)
			rect(data[i,6], 28, data[i,7], 32, , col = "blue")
			rect(data[i,4] + gap, 28, data[i,5] + gap, 32, col = "blue")
			lines(c(data[i,7], data[i,4] + gap), c(30,30))
			lines(c(data[i,6], data[i,10]), c(28, 12), lty = 2)
                        lines(c(data[i,7], data[i,11]), c(28, 12), lty = 2)
                        lines(c(data[i,4] + gap, data[i,8]), c(28, 12), lty = 2)
                        lines(c(data[i,5] + gap, data[i,9]), c(28, 12), lty = 2)
			text(c(data[i,4] + gap, data[i,5] + gap, data[i,6], data[i,7]), c(33, 33, 33, 33), c(data[i,4], data[i,5], data[i,6], data[i,7]))                        
			text(c(data[i,8], data[i,9], data[i,10], data[i,11]), c(7, 7, 7, 7), c(data[i,8], data[i,9], data[i,10], data[i,11]))
	}
	text(max(data[i,4:11])/2, 35, data[i,1], cex = 2, font = 2)
	text(max(data[i,4:11])/2, 15, data[i,2], cex = 2, font = 2)
	rect(min(data[i,8:11]), 8, max(data[i,8:11]), 12, col = "blue")
	dev.off()
}

data <- data[, -3]
colnames(data) <- c("Transcript1", "Transcript2", "T1_Start1", "T1_End1", "T1_Start2", "T1_End2", "T2_Start1", "T2_End1", "T2_Start2", "T2_End2")
write.table(data, paste(opt$o, "/as_summary.xls", sep = ""), sep = "\t", row.names = F, quote = F)
