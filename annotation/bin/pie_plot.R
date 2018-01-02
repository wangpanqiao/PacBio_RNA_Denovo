# usage function
usage <- function() {
        print("-------------------------------------------------------------------------------")
        print("Usage: Rscript go_terms_stat.r stat.xls out_prefix database graph_lable")
        print("1) stat.xls: the annotation stat results (eg. Evalue & Identity & Species)")
        print("2) out_prefix: the png & pdf files prefix")
        print("3) database: the database of stat results")
        print("4) graph_lable: the graph title of stat results (eg. Evalue & Identity & Species)")
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
database <- args[3]
title <- args[4]

.libPaths("/home/mengfei/workdir/R_install_packges/")

df <- read.table(file=infiles, header = FALSE, quote = "\"'", dec = ".", sep = "\t", row.names = 1, comment.char = "#")

png (filename = paste (outfile, ".png", sep = ""), width=1000, height = 750, res = 30, units = "px")
pie (df[,1], main = paste ("Blast", title, "Distribution", "(", database, ")", sep = " "), col = rainbow(length(df[,1])), labels = df[,2], clockwise = TRUE, cex.main = 5.0, cex = 3.0)
legend ("right", row.names(df), cex = 3.0, fill = rainbow(length(df[,1])))
#print(p)
dev.off ()

pdf (file = paste (outfile, ".pdf", sep = ""), width = 12, height = 8)
pie (df[,1], main = paste ("Blast", title, "Distribution", "(", database, ")", sep = " "), col = rainbow(length(df[,1])), labels = df[,2], clockwise = TRUE, cex.main = 1.5, cex = 0.8)
legend ("right", row.names(df), cex = 1.0, fill = rainbow(length(df[,1])))
#print(p)
dev.off ()
