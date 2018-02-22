# Copyright ---------------------------------------------------------------
# 2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# signal to noise scoring function

# Input -------------------------------------------------------------------
# Gene expression data and genesets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------


# Library imports ---------------------------------------------------------

# Determine signal to noise scoring for each line in the data table.  dataA is
# the baseline(negative if higher) and dataB is the test(positive if higher)

#' Title
#'
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
s2n.score <- function(a, b) {
  mu1 <- apply(a, 1, mean)
  mu2 <- apply(b, 1, mean)
  sigma1 <- apply(a, 1, sd)
  sigma2 <- apply(b, 1, sd)
  (mu2 - mu1) / (sigma1 + sigma2)
}
