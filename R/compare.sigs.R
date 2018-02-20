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
# library(parallel)
# library(reshape2)

#' Compares sparse matrixes of gene signatures
#' 
#' @param sig.matrix Signature sparse matrix
#' @param test.matrix Test signatures sparse matrix
#'   
#' @return A data table of the comparisons of sig.matrix and test.matrix
#'   including odds ratios, p.values and FDR adjusted q.values
#' @examples
#' compare.sigs(sig.matrix, test.matrix)
compare.sigs <- function(sig.matrix, test.matrix) {
  # make the comparison
  system.time(results.vec <- sig.matrix %*% t(test.matrix))
  # get the other parts we need
  genes.per.test <- rowSums(test.matrix)
  genes.per.sig <- rowSums(sig.matrix)
  total.genes <- ncol(sig.matrix)
  gc()
  library(parallel)
  system.time(p.values <- mcmapply(function(i, j) {
    phyper(results.vec[i, j], genes.per.sig[i],
           total.genes - genes.per.sig[i], genes.per.test[j], F)
  },
  i = rep(1:nrow(sig.matrix), times = length(genes.per.test)),
  j = rep(1:length(genes.per.test), each = nrow(sig.matrix)), mc.cores = 10))
  
  
  dim(p.values) <- c(nrow(sig.matrix), length(genes.per.test))
  rownames(p.values) <- row.names(sig.matrix)
  colnames(p.values) <- names(genes.per.test)
  odds.ratio <- total.genes * results.vec / genes.per.sig ^ 2
  
  
  library(reshape2)
  output <- reshape2::melt(odds.ratio)
  output2 <- reshape2::melt(p.values)
  output$p <- output2$value
  colnames(output) <- c("sig", "test.sig", "odds.ratio", "p.value")
  output$q.value <- p.adjust(output$p.value, "BH")
  
  output
}