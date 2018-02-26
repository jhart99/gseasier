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
#'
#' @return
#' @export
#'
#' @examples
FES <- function(sn.table, geneset, hits.only=TRUE) {
  genes <- rownames(sn.table)
  overlap <- intersect(genes, geneset)
  hit <- genes %in% overlap
  # preallocate
  lower <- vector("numeric", nrow(sn.table))
  higher <- vector("numeric", nrow(sn.table))
  hit.count <- 0
  miss.count <- 0
  for (index in 1:nrow(sn.table)) {
    if (hit[index]) {
      hit.count <- hit.count + 1
      lower[index] <- phyper(hit.count, length(overlap),
                             length(genes) - length(overlap), index)
      higher[index] <- phyper(length(overlap) - hit.count, length(overlap),
                              length(genes) - length(overlap),
                              length(genes) - index)
    } else {
      miss.count <- miss.count + 1
      if (hits.only) {
        lower[index] <- 1
        higher[index] <- 1
      } else {
        lower[index] <- phyper(hit.count, length(overlap),
                               length(genes) - length(overlap), index)
        higher[index] <- phyper(length(overlap) - hit.count,
                                length(overlap),
                                length(genes) - length(overlap),
                                length(genes) - index)
      }
    }
  }
  lower[lower > 1] <- 1
  higher[higher > 1] <- 1
  data.frame(sn.table, p.lower = lower, p.higher = higher, hit = hit)
}
