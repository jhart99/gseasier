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
z.score <- function(gct, cls) {
  exp.matrix <- as.matrix(gct[, c(-1, -2)])
  exp.matrix <- log2(exp.matrix)
  mean.int <- apply(exp.matrix, 2, mean)
  sd.int <- apply(exp.matrix, 2, sd)
  z.score.mat <- foreach(i=1:ncol(exp.matrix), .combine=cbind) %do% {
    exp.matrix[, i] <- (exp.matrix[, i] - mean.int[i])/sd.int[i]
  }
  colnames(z.score.mat) <- colnames(exp.matrix)
  data.frame(gct[, c(1,2)], z.score.mat)
  
}