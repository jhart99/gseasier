# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Enrichment set plot similar to the one provided by the GSEA program 

# Input -------------------------------------------------------------------


# Methods -------------------------------------------------

# Outputs -------------------------------------------------

# Library imports ---------------------------------------------------------

#' Enrichment set plot
#'
#' @param gmt A gmt formatted gene expression data set
#' @param cls A cls formatted classifier
#' @param comparison Which levels of cls should be compared defaults to levels 1 and 2
#' @param geneset Either a vector of gene identifiers or a geneset name in combination with msigdb parameter
#' @param sn.table (optional) A signal to noise named vector can be used in place of gmt, cls and comparison
#' @param msigdb (optional) A database of signatures
#'
#' @return
#' @export
#'
#' @examples
ES.plot <- function(gmt=NULL, cls=NULL, comparison=NULL,
                    geneset=NULL, sn.table=NULL, msigdb=NULL) {
  library(grid)
  library(ggplot2)
  library(gtable)
  library(gridExtra)
  source("R/ES.R")

  if (is.null(gmt)) {
    if (is.null(sn.table)) {
      stop("Either gmt and cls or sn.table must be used to generate a plot")
    }
  } else {
    if (length(comparison) == 2 | is.null(comparison)) {
      sn.table <- s2n.rank(gmt, cls, comparison)
    } else {
      stop("Invalid comparison")
    }
  }

  if (length(geneset) == 1) {
    if (is.null(msigdb)) {
      stop("msigdb missing")
    }
    geneset <- msigdb$gene[msigdb$sig == geneset]
  }

  plot.data <- ES(sn.table, geneset, hits.only = F)

  plot.data$percentile <- cut(plot.data$sn,
                              quantile(plot.data$sn, probs = seq(0, 1, 1 / 9)))

  grid.newpage()
  sn.plot <- ggplot(plot.data, aes(x = rank, y = sn)) +
    geom_point() + theme_bw() +
    theme(plot.margin = unit(c(0, 1, 0, 0), "lines"),
          panel.border = element_rect(size = 0),
          panel.grid = element_blank())
  ES.plot <- ggplot(plot.data, aes(x = rank, y = ES)) +
    geom_line(color = "blue") + theme_bw() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0, 1, 0, 0), "lines"),
          panel.border = element_rect(size = 0),
          panel.grid = element_blank())
  ticks <- ggplot(plot.data, aes(x = rank,
                                 xend = rank,
                                 y = ifelse(hit, 0.1, 0),
                                 yend = 0)) +
    geom_segment() + theme_bw() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 1, 0, 1.1), "lines"),
          panel.border = element_rect(size = 0),
          panel.grid = element_blank())
  gradient <- ggplot(plot.data, aes(x = rank)) + theme_bw() +
    geom_tile(aes(y = 0, fill = as.numeric(percentile))) +
    scale_fill_gradient2(midpoint = 5, low = "blue", high = "red") +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 1, 0, 1.1), "lines"),
          panel.border = element_rect(size = 0),
          panel.grid = element_blank())

    # hacky things to modify the facets to different heights

  grid.arrange(ES.plot, ticks, gradient, sn.plot,
               ncol = 1, nrow = 4, heights = c(5, 1, 1, 3))
}
