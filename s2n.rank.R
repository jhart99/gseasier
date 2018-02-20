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
#' @param cls 
#' @param contrast 
#'
#' @return
#' @export
#'
#' @examples
s2n.rank <- function(gct, cls, contrast=NULL) {
  if (is.null(contrast)){
    contrast=list(levels(cls)[1], levels(cls)[2])
  }
  sn.vector <- s2n.score(gct[, which(cls %in% contrast[[1]])+2], gct[, which(cls %in% contrast[[2]])+2])
  sn.vector <- data.frame(sn=sn.vector)
  rownames(sn.vector) <- gct[, 1]
  sn.vector <- na.omit(sn.vector)
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector <- sn.vector[order(-sn.vector$sn), ]
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector
}