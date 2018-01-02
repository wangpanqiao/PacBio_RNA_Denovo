#!/share/public/software/R-3.3.3/bin/Rscript
# coding=utf-8
# Writer:         chaijingchao
# Program Date:   2016.06.30

library(ggplot2)
args <- commandArgs(TRUE)
length_data <-read.table(args[1])
names(length_data) <- c("cDNA_size", "read_id", "read_length")
p <- ggplot(length_data, aes(read_length, colour = cDNA_size, fill = cDNA_size)) + geom_density(alpha = 0.1) +  xlab("Sequence length(bp)") + xlim(0, 10000)
#p <- p + scale_color_brewer("type")
p <- p + theme(panel.border = element_rect(fill = NA,colour = "black"))
p <- p + theme(
panel.background = element_rect(fill = "transparent",colour =NA),
panel.grid.minor = element_blank(),
panel.grid.major = element_blank(),
plot.background = element_rect(fill  = "transparent",colour =NA)
)
ggsave(plot = p, file = paste(args[2], ".pdf", sep=""), width = 6, height = 4)
ggsave(plot = p, file = paste(args[2], ".png", sep=""), width = 6, height = 4, type = 'cairo-png')
