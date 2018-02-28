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
FES.summary <- function(sn.table, geneset) {
  fes.data <- FES(sn.table, geneset)
  high.point <- which.min(fes.data$p.higher)
  low.point <- which.min(fes.data$p.lower)
  with(fes.data, data.frame(overlap = sum(hit),
                            max.lo = low.point,
                            sn.lo = sn[low.point],
                            p.lo = p.lower[low.point],
                            max.hi = high.point,
                            sn.hi = sn[high.point],
                            p.hi = p.higher[high.point]))
}
