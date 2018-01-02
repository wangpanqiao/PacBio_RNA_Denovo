# usage function
usage <- function() {
        print("-------------------------------------------------------------------------------")
        print("Usage: Rscript go_terms_stat.r stat.xls out_prefix database")
        print("1) stat.xls: the annotation stat results (eg. Evalue & Identity)")
        print("2) out_prefix: the png & pdf files prefix")
        print("3) database: the database of stat results")
        print("-------------------------------------------------------------------------------")
}

# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) < 3 ) {
        print(args)
        usage()
        stop("the length of args < 2")
}

infiles <- args[1]
outfile <- args[2]
database <- args[3]

.libPaths("/home/mengfei/workdir/R_install_packges/")

df <- read.table(file=infiles, header = FALSE, quote = "\"'", dec = ".", sep = "\t", row.names = 1, comment.char = "#")

png (filename = paste (outfile, ".png", sep = ""), width=1400, height = 800, res = 30, units = "px")
pie (df[,1], main = paste ("Blast Species Distribution", "(", database, ")", sep = " "), col = rainbow(length(df[,1])), labels = df[,2], clockwise = TRUE, cex.main = 5.0, cex = 3.0)
legend ("right", row.names(df), cex = 3.0, fill = rainbow(length(df[,1])))
#print(p)
dev.off ()

pdf (file = paste (outfile, ".pdf", sep = ""), width = 14, height = 9)
pie (df[,1], main = paste ("Blast Species Distribution", "(", database, ")", sep = " "), col = rainbow(length(df[,1])), labels = df[,2], clockwise = TRUE, cex.main = 1.5, cex = 0.8)
legend ("right", row.names(df), cex = 1.0, fill = rainbow(length(df[,1])))
#print(p)
dev.off ()
