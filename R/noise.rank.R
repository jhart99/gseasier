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
noise.rank <- function(gct, cls, contrast=NULL) {
  if (is.null(contrast)) {
    contrast <- list(levels(cls)[1])
  }
  sd.vector <- apply(gct[, which(cls %in% contrast[[1]]) + 2], 1, sd)
  sn.vector <- data.frame(sn = sd.vector)
  rownames(sn.vector) <- gct[, 1]
  sn.vector <- na.omit(sn.vector)
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector <- sn.vector[order(-sn.vector$sn), ]
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector
}
