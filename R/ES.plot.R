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
  
  # hacky things to modify the facets to different heights, remove the axis that
  # aren't needed, etc. this part is fragile and might break if ggplot2 is a
  # different version.  These are magic numbers here and I will try to explain
  # what they are, how to get them and how to fix this if ggplot updates and
  # blows this part up.
  
  # in order to do advanced things to the plot we need to make it into a grid 
  # object so we can adjust the parts. First 15 and 16 are targetting the Y axes
  # for panel 2 and 3 and resetting them to null.  Second, 7 and 22 are 
  # adjusting the heights of the panel 1 and panel 4.  The other parts are a 
  # height of 1 so top is 3 times and bottom is 2 times of the middle parts. 
  # Finally add a second Axis label.  This is done here because if we put this
  # in ggplot the panels don't align right.  the 7 and 22 refer to the positions
  # of panel 1 and panel 4 in the table.
  
  # how to get to the numbers?  well set a breakpoint in this function. after 
  # the following line.  look at gt and you will see all the pieces of the plot.
  # Beware, the order doesn't correspond to the numbers you want.  Numbers for
  # the second axis are easy just look for these.
  
  # 2   1 ( 7- 7, 4- 4)   panel-1-1               gTree[panel-1.gTree.6101]
  # 3   1 (12-12, 4- 4)   panel-1-2               gTree[panel-2.gTree.6116]
  # 4   1 (17-17, 4- 4)   panel-1-3               gTree[panel-3.gTree.6131]
  # 5   1 (22-22, 4- 4)   panel-1-4               gTree[panel-4.gTree.6146]
  
  # Those tell you to use 7 and 22 for the text height and you need to use
  # something less than 4 for the width.  3 is the normal axis, so go with 2.
  
  # same for the axis to delete
  
  # 14  3 ( 7- 7, 3- 3)  axis-l-1-1    absoluteGrob[GRID.absoluteGrob.5436]
  # 15  3 (12-12, 3- 3)  axis-l-2-1    absoluteGrob[GRID.absoluteGrob.5443]
  # 16  3 (17-17, 3- 3)  axis-l-3-1    absoluteGrob[GRID.absoluteGrob.5450]
  # 17  3 (22-22, 3- 3)  axis-l-4-1    absoluteGrob[GRID.absoluteGrob.5457]
  
  # get those and set them to zeroGrob()
  
  # finally you need to adjust the heights, but these don't correspond to the 
  # same order.  Basically you will have to do this by inspecting the height 
  # object. gt$heights These are ugly but you are looking for the lines that
  # have a height of 1null those are the panel heights.  They are in order from
  # top to bottom.
  
  # [1] 0lines              0cm                 0cm                 0cm                 0cm                 0cm                 1null               0cm                 0lines              0cm                
  # [11] 0cm                 1null               0cm                 0lines              0cm                 0cm                 1null               0cm                 0lines              0cm                
  # [21] 0cm                 1null               0.412097602739726cm 1grobheight         0cm                 0lines  
  
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