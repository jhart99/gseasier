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
#'
#' @return
#' @export
#'
#' @examples
leading.edge <- function(sn.table, geneset) {
  scores <- ES(sn.table, geneset)
  peak <- which.max(abs(scores$ES))
  if (scores$ES[peak] > 0) {
    return(rownames(scores)[1:peak])
  } else {
    return(rownames(scores)[peak:length(scores)])
  }
}
