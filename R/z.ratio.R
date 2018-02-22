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
z.ratio <- function(a, b) {
  mu1 <- apply(a, 1, mean)
  mu2 <- apply(b, 1, mean)
  sigma1 <- apply(a, 1, sd)
  sigma2 <- apply(b, 1, sd)
  n1 <- ncol(a)
  n2 <- ncol(b)
  (mu1 - mu2) / sqrt( (sigma1 ^ 2 / n1) + (sigma2 ^ 2 / n2))
}
