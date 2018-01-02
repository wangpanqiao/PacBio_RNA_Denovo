plot_MA = function(logCounts, logFoldChange, FDR, FDR_thresh, log2FC_thresh, xlab="logCounts", ylab="logFC", title="MA plot", pch=20) {

    plot(logCounts, logFoldChange, col=ifelse(FDR > FDR_thresh, "blue",ifelse(logFoldChange >= log2FC_thresh,"red",ifelse(logFoldChange <= -log2FC_thresh, "green","blue"))), xlab=xlab, ylab=ylab, main=title, pch=pch);;

}


plot_Volcano = function(logFoldChange, FDR, FDR_thresh, log2FC_thresh, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20) {

   fdr_adj = ifelse(FDR<=1e-50,1e-50,FDR);
   plot(logFoldChange, -1*log10(fdr_adj), col=ifelse(fdr_adj>FDR_thresh, "blue",ifelse(logFoldChange >= log2FC_thresh,"red",ifelse(logFoldChange <= -log2FC_thresh, "green","blue"))), xlab=xlab, ylab=ylab, main=title, pch=pch);

}


plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, FDR_thresh, log2FC_thresh, xlab="logCounts", ylab="logFC", title="MA plot") {

    def.par = par(no.readonly = TRUE) # save default, for resetting...

    gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
    layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 

    plot_MA(logCounts, logFoldChange, FDR, FDR_thresh, log2FC_thresh);
    plot_Volcano(logFoldChange, FDR, FDR_thresh, log2FC_thresh);

    # draw again, but use a smaller dot for data points
    plot_MA(logCounts, logFoldChange, FDR, FDR_thresh, log2FC_thresh, pch='.');
    plot_Volcano(logFoldChange, FDR, FDR_thresh, log2FC_thresh, pch='.');
    

    par(def.par)   
        
    
}