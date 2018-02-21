# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Converts the data.frame used by msigdb to a sparse matrix

# Input -------------------------------------------------------------------
# 

# Methods -------------------------------------------------
# convert.sigs.to.matrix

# Outputs -------------------------------------------------
# a sparse matrix of genes by sets

# Library imports ---------------------------------------------------------

#' Converts signatures into a sparse matrix
#' 
#' @param signatures Signatures table from loadSig
#' @param genes A list of all the genes that will be used in the comparison
#'
#' @return A sparse matrix of signatures by gene
#' @examples
#' convert.sigs.to.matrix(signatures, genes)
convert.sigs.to.matrix <- function(signatures, genes) {
  # This could be rewritten using some sort of sparse matrix which will
  # save on memory, but in general I haven't found that the performance tradeoff
  # is worth it.  Using Matrix reduces memory usage by ~10 fold but increases
  # CPU time by ~10x or more.

  # gene names will be used as column names so we need to make sure they are R
  # safe
  signatures$xgene <- make.names(signatures$gene)
  signatures$xsig <- make.names(signatures$sig)
  total.genes <- length(genes)
  sig.labels <- unique(signatures$xsig)
  sig.matrix <- matrix(0L, nrow = length(sig.labels), ncol = total.genes)
  rownames(sig.matrix) <- sig.labels
  colnames(sig.matrix) <- genes
  system.time(for (i in 1:nrow(signatures)) {
    sig.matrix[signatures$xsig[i], signatures$xgene[i]] <- 1L
  })
  sig.matrix
}