# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A demonstration of the gene signature overlap routines.

# Input -------------------------------------------------------------------
# msigdb formatted gmt files or yuor own gmt formatted gene sets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------

# a table of comparisons between two sets of gene sets with odds ratios and
# p.values

# Library imports ---------------------------------------------------------

# acquire the msigdb files from Broad directly at
# http://software.broadinstitute.org/gsea/downloads.jsp

# this part will work with either set of genes.  Meaning you can use entrez ids
# or gene names as you prefer.


#' Title
#'
#' @param sn.table
#' @param geneset
#' @param hits.only
#' @param weight
#'
#' @return
#' @export
#'
#' @examples
ES <- function(sn.table, geneset, hits.only=T, weight=1) {
  genes <- rownames(sn.table)
  overlap <- intersect(genes, geneset)
  hit <- genes %in% overlap

  hit.vals <- abs(sn.table$sn  ** weight) * hit

  hit.scale <- 1 / sum(hit.vals)
  miss.scale <- 1 / (length(genes) - length(overlap))

  # preallocate
  out <- vector("numeric", nrow(sn.table))

  current.score <- 0
  for (index in 1:nrow(sn.table)) {
    score <- ifelse(hit[index],
                    hit.vals[index] * hit.scale, -miss.scale) + current.score
    out[index] <- score
    current.score <- score
  }
  data.frame(sn.table, ES = out, hit = hit)
}
