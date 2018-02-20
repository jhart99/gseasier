# Copyright ---------------------------------------------------------------
# 2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Filter a gct object for genes with low expression

# Input -------------------------------------------------------------------
# Gene expression data and genesets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------


# Library imports ---------------------------------------------------------



#' Title
#'
#' @param gct 
#' @param cls 
#' @param min.cpm 
#' @param contrast 
#'
#' @return
#' @export
#'
#' @examples
filter.cpm <- function(gct, cls, min.cpm=1, contrast=NULL) {
  if (is.null(contrast)){
    contrast=list(levels(cls)[1], levels(cls)[2])
  }
  len1 <- sum(cls %in% contrast[[1]])
  len2 <- sum(cls %in% contrast[[2]])
  selection1 <- rowSums(gct[, which(cls %in% contrast[[1]])+2]>min.cpm)
  selection2 <- rowSums(gct[, which(cls %in% contrast[[2]])+2]>min.cpm)
  gct[(selection1==len1) & (selection2==len2),]
}