# coding=utf-8
# Writer:         chaijingchao
# Program Date:   2016.07.30

library(ggplot2)
args <- commandArgs(TRUE)
data <- read.table("lnc_length.xls")
names(data) <- c("cDNA","gene_id","length")
p <- ggplot(data,aes(length)) + geom_histogram(binwidth=500, colour = "black", fill = "steelblue")
#p <- p + ggtitle("Length Distribution") 
p <- p + theme(panel.border = element_rect(fill = NA,colour = "black"))
p <- p + theme(
panel.background = element_rect(fill = "transparent",colour =NA),
panel.grid.minor = element_blank(),
panel.grid.major = element_blank(),
plot.background = element_rect(fill  = "transparent",colour =NA)
)
ggsave(plot = p, file = paste(args[2], ".pdf", sep=""), width = 6, height = 4)
ggsave(plot = p, file = paste(args[2], ".png", sep=""), width = 6, height = 4, type = 'cairo-png')

