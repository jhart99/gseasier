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
#'
#' @return
#' @export
#'
#' @examples
z.scale <- function(gct, cls) {
  exp.matrix <- as.matrix(gct[, c(-1, -2)])
  exp.matrix <- log2(exp.matrix)
  # row medians
  median.vector <- apply(exp.matrix, 1, median)
  # sd by sample group and z-scale
  scaled.matrix <- foreach(sample.group=iter(unique(cls)), .combine=cbind) %do% {
    sd.vector <- apply(exp.matrix[, which(cls == sample.group)], 1, sd)
    (exp.matrix[, which(cls == sample.group)] - median.vector)/sd.vector
    #     sd.vector <- data.frame(sd.vector)
    #     colnames(sd.vector) <- sample.group
    #     sd.vector
  }
  data.frame(gct[, c(1,2)], scaled.matrix)
  
}