# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A gene set overlap functions

# Input -------------------------------------------------------------------
# convert.sigs.to.matrix:
# msigdb formatted gmt files or yuor own gmt formatted gene sets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------
# a sparse matrix of genes by sets


# Library imports ---------------------------------------------------------

#' Title
#'
#' @param gsea
#'
#' @return
#' @export
#'
#' @examples
gsea.pkg.env <- new.env(parent = emptyenv())
set.gsea <- function(gsea) {
  if (file.exists(gsea)) {

    gsea.pkg.env$gsea.exec <- gsea
  } else {
    stop("invalid gsea executable")
  }
}
