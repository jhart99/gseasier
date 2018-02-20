# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A demonstration of the gene signature overlap routines. 

# Input -------------------------------------------------------------------
# msigdb formatted gmt files or yuor own gmt formatted gene sets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------

# a table of comparisons between two sets of gene sets with odds ratios and
# p.values

# Library imports ---------------------------------------------------------

# acquire the msigdb files from Broad directly at
# http://software.broadinstitute.org/gsea/downloads.jsp

# this part will work with either set of genes.  Meaning you can use entrez ids
# or gene names as you prefer.


#' Title
#'
#' @param gmt 
#' @param cls 
#' @param comparison 
#' @param geneset 
#'
#' @return
#' @export
#'
#' @examples
ES.plot <- function(gmt, cls, comparison, geneset) {
  library(ggplot2)
  library(grid)
  library(gtable)
  
  source('R/s2n.rank.R')
  source('R/ES.R')
  sn.table <- s2n.rank(gmt, cls, comparison)
  plot.data <- ES(sn.table, geneset, hits.only=F) 
  plot.data$percentile <- as.numeric(cut(plot.data$sn, quantile(plot.data$sn, probs=seq(0, 1, 1/9)), include.lowest=T))
  plot.data$hit <- ifelse(plot.data$hit, 1, 0)
  
  peak.pos <- plot.data$rank[which.max(abs(plot.data$ES))]
  
  plot.data2 <- melt(plot.data, id.vars="rank")
  
  # fix the facet order
  plot.data2$variable <- factor(plot.data2$variable, levels=c("ES", "hit", "percentile", "sn"))
  
  p <- ggplot(plot.data2, aes(x=rank)) + 
    facet_wrap(~ variable, ncol=1, scales="free_y") +
    geom_line(data=plot.data2[plot.data2$variable == "ES", ], aes(y=value), color="darkgreen") +
    
    geom_point(data=plot.data2[plot.data2$variable == "sn", ], aes(y=value)) +
    geom_segment(data=plot.data2[plot.data2$variable == "hit", ], aes(x=rank, xend=rank, y=value, yend=0)) +
    geom_tile(data=plot.data2[plot.data2$variable == "percentile", ], aes(y=1, fill=value)) + scale_fill_gradient2(midpoint=5, low="blue", high="red") +
    #geom_rect(aes(xmin=index.min, xmax=index.max, color=index), ymin=0, ymax=1, data=percentile.data) +
    #geom_segment(subset=.(variable == "percentile"), aes(y=1, yend=0, x=rank, xend=rank, color=value)) + scale_color_gradient2(midpoint=5, low="blue", high="red") +
    geom_hline(yintercept=0, linetype=1) +
    geom_vline(xintercept=peak.pos) +
    theme_bw() +
    theme(strip.text=element_blank(), 
          legend.position="none", 
          strip.background = element_blank(),
          plot.margin = unit(c(0,0,0,0) , units = "lines"),
          panel.spacing= unit(c(0), units="lines"),
          panel.border=element_blank()) +
    ylab("")
  
  # hacky things to modify the facets to different heights, remove the axis that aren't needed, etc.
  # this part is fragile and might break if ggplot2 is a different version.
  gt <- ggplot_gtable(ggplot_build(p))
  gt$grobs[[15]] <- zeroGrob()
  gt$grobs[[16]] <- zeroGrob()
  gt$heights[[7]] <- unit(3, "null") 
  gt$heights[[22]] <- unit(2, "null")
  gt <- gtable_add_grob(gt, rectGrob(gp=gpar(fill="white", col="white")), 7, 7, 22, name="blank1")
  gt <- gtable_add_grob(gt, textGrob("Enrichment Score", rot=90, gp=gpar(fontsize=11)), 7, 2, name="blank1", clip="off")
  gt <- gtable_add_grob(gt, textGrob("Signal to Noise", rot=90, gp=gpar(fontsize=11)), 22, 2, name="blank1", clip="off")
  grid.draw(gt) 
}


#' Title
#'
#' @param sn.table 
#' @param geneset 
#'
#' @return
#' @export
#'
#' @examples
ES.plot2 <- function(sn.table, geneset) {
  library(grid)
  library(ggplot2)
  library(gtable)
  library(gridExtra)
  source('R/ES.R')
  plot.data <- ES(sn.table, geneset, hits.only=F) 
  
  plot.data$percentile <- cut(plot.data$sn, quantile(plot.data$sn, probs=seq(0, 1, 0.11111)))
  #ggplot(plot.data, aes(x=rank)) + geom_line(aes(y=-log10(p.lower)), color="blue") + geom_line(aes(y=-log10(p.higher)), color="red")
  
  grid.newpage()
  sn.plot <- ggplot(plot.data, aes(x=rank, y=sn)) + geom_point() + theme_bw() +
    theme(plot.margin=unit(c(0, 1, 0, 0), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  ES.plot <- ggplot(plot.data, aes(x=rank, )) + geom_line(aes(y=-log10(p.lower)), color="blue") +
    geom_line(aes(y=-log10(p.higher)), color="red") +
    theme_bw() +
    theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(plot.margin=unit(c(0, 1, 0, 0), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  ticks <- ggplot(plot.data, aes(x=rank, xend=rank, y=ifelse(hit, 0.1, 0), yend=0)) +geom_segment() + theme_bw() +
    theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.line=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(plot.margin=unit(c(0, 1, 0, 1.1), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  gradient <- ggplot(plot.data, aes(x=rank))  +theme_bw() +
    geom_tile(aes(y=0, fill=as.numeric(percentile))) + scale_fill_gradient2(midpoint=5, low="blue", high="red") + theme(legend.position="none") +
    theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.line=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(plot.margin=unit(c(0, 1, 0, 1.1), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  
  # hacky things to modify the facets to different heights
  
  
  grid.arrange(ES.plot, ticks, gradient, sn.plot, ncol=1, nrow=4, heights=c(5, 1, 1, 3))
  
  
}