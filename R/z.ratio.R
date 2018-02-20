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
#' @param dataA 
#' @param dataB 
#'
#' @return
#' @export
#'
#' @examples
z.ratio <- function(dataA, dataB) {
  mu1 <- apply(dataA, 1, mean)
  mu2 <- apply(dataB, 1, mean)
  sigma1 <- apply(dataA, 1, sd)
  sigma2 <- apply(dataB, 1, sd)
  n1 <- ncol(dataA)
  n2 <- ncol(dataB)
  (mu1-mu2)/sqrt((sigma1^2/n1)+(sigma2^2/n2))
}