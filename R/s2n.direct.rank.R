# Copyright ---------------------------------------------------------------
# 2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Fisher enrichment and accessory functions

# Input -------------------------------------------------------------------
# Gene expression data and genesets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------


# Library imports ---------------------------------------------------------

#' Title
#'
#' @param gct 
#' @param high.samples 
#' @param low.samples 
#'
#' @return
#' @export
#'
#' @examples
s2n.direct.rank <- function(gct, high.samples, low.samples) {
  sn.vector <- s2n.score(gct[, which(cls %in% contrast[[1]])+2], gct[, which(cls %in% contrast[[2]])+2])
  sn.vector <- data.frame(sn=sn.vector)
  rownames(sn.vector) <- gct[, 1]
  sn.vector <- na.omit(sn.vector)
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector <- sn.vector[order(-sn.vector$sn), ]
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector
}