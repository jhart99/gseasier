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
#' @param files
#'
#' @return
#' @export
#'
#' @examples
read.gsea.results <- function(files) {
  results <- do.call(rbind, lapply(files, function(file) {
    html.in  <- XML::htmlParse(file)
    results.table <- XML::readHTMLTable(html.in, header = F,
                                        as.data.frame = F, which = 1)
    results.table <- data.frame(lapply(results.table, type.convert, as.is = T))
    results.table <- results.table[, c(-1, -3)]
    colnames(results.table) <- c("sig", "size", "ES",
                                 "NES", "p", "q",
                                 "fwer", "rank", "leader")
    results.table
  }))
  results[order(results$ES), ]
}
