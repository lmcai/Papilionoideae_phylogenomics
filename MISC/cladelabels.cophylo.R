# New function to plot clade labels from a cophylo object
# Author: Domingos Cardoso

clabelabels.cophylo <- function(cophylo,
                                label = NULL,
                                taxa = NULL, 
                                which= c("left","right"),
                                col="black",
                                cex=1,
                                offset=NULL) {
  
  lastPP <- get("last_plot.cophylo", envir=.PlotPhyloEnv)
  
  if(is.null(offset)) offset <- 0.02
  
  if (which[1]=="left") {
    pos_label_left <- which(cophylo$trees[[1]]$tip.label %in% taxa)
    lx <- lastPP$left$xx[pos_label_left]
    ly <- lastPP$left$yy[pos_label_left]
    text(lx+offset, ly, label, cex=cex, col=col, adj=c(0, 0.5))
    
  } else if(which[1]=="right") {
    pos_label_right <- which(cophylo$trees[[2]]$tip.label %in% taxa)
    rx <- lastPP$right$xx[pos_label_right]
    ry <- lastPP$right$yy[pos_label_right]
    text(rx-offset, ry, label, cex=cex, col=col, adj=c(1, 0.5))
  }
  
}

