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
#' @param doc.string
#' @param output.directory
#'
#' @return
#' @export
#'
#' @examples
find.gsea.results <- function(doc.string, output.directory="gsea") {
  report.dir <- list.files(
    paste(getwd(), output.directory, doc.string, sep = "/"),
    "gsea_report.*html", full.names = T)
  report.dir
}
