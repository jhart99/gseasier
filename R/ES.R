# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Calculates the enrichment score for a geneset vs a ranking

#' Enrichment score
#'
#' @param sn.table an ordered, valued ranking of genes
#' @param geneset a set of genes to test
#' @param hits.only Emit values only for genes in the geneset(True) or all values(False)
#' @param weight The weighting parameter used in GSEA v2 and up.  For V1 this value is set to 0.
#'
#' @return A data.frame of the enrichment score for each gene in the geneset or ranking.
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
  if (hits.only) {
    data.frame(sn.table[hit, ], ES = out[hit], hit = hit[hit])
  } else {
    data.frame(sn.table, ES = out, hit = hit)
  }
}
