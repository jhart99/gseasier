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
#' @param genesets
#'
#' @return
#' @export
#'
#' @examples
FES.table <- function(sn.table, genesets) {
  # this factor is to drop unused levels from the sig list which causes an issue
  # with ddply
  genesets$sig <- factor(genesets$sig)
  out.table <- ddply(genesets,
                     .(sig),
                     .parallel = T,
                     .fun = function(df) FES.summary(sn.table, df$gene))
  out.table$q.lo <- p.adjust(out.table$p.lo, method = "BH")
  out.table$q.hi <- p.adjust(out.table$p.hi, method = "BH")
  out.table[, c(1:5, 9, 6:8, 10)]
}
