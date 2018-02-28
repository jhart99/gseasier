# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A gene set overlap function

#' Title
#'
#' @param sn
#'
#' @return
#' @export
#'
#' @examples
jitter.sn <- function(sn) {
  sn.mean <- mean(sn$sn)
  sn.sd <- sd(sn$sn)
  sn$sn <- sn$sn + rnorm(sn, sn.mean, sn.sd)
  sn <- sn[order(-sn$sn), ]
  sn$rank <- 1:nrow(sn)
  sn
}
