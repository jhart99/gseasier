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
#' @param sigs
#'
#' @return
#' @export
#'
#' @examples
leading.edge.heatmap <- function(sn.table, sigs) {
  sigs.to.test <- unique(sigs$sig)
  leading.genes <- lapply(sigs.to.test, function(x) {
    leading.edge(sn.table, as.character(sigs$gene[sigs$sig == x]))
  })
  leading.gene.table <- data.frame(ragged.list.melt(leading.genes))
  leading.gene.table$i <- sigs.to.test[leading.gene.table$i]
  colnames(leading.gene.table) <- c("sig", "gene")

  # you need to know the set of all genes present in the signatures and your test
  # set.  We need to know this because in the comparison stage the matrixes need
  # the same dimensions and order.
  genes <- make.names(unique(leading.gene.table$gene))
  # convert the gene signature table into a sparse matrix
  test.matrix <- convert.sigs.to.matrix(leading.gene.table, genes)
  test.compare <- compare.sigs(test.matrix, test.matrix)

  library(gplots)
  library(reshape2)
  test.comp.matrix <- acast(test.compare, sig ~ test.sig,
                            value.var = "q.value")
  test.comp.matrix[test.comp.matrix == 0] <- 1e-100
  test.comp.matrix <- -log2(test.comp.matrix)
  heatmap.2(test.comp.matrix, trace = "none", scale = "none",
            margins = c(15,15))
}
