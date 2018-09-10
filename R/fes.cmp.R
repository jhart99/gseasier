# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A demonstration of the gene signature overlap routines.




#' FES
#'
#' FES(Fast Enrichment Score or Fisher Enrichment Score) uses a different method
#' of calculating the enrichment score as compared with GSEA.  In FES, the
#' hypergeometric distribution is used to determine the p-value for the
#' partitioning of members of the geneset within the ordered valued gene
#' ranking.  Two different p-values are calculated from the opposite ends of the
#' ranking to measure enrichment at either the high or low ends of the spectrum.
#' The use of the hypergeometric distribution rather than the complicated method
#' used in GSEA means this method is very fast and yields results which are
#' similar to GSEA itself.  In comparison to GSEA, FES is more sensitive to
#' concentration of geneset members in the middle of the ranking as happens with
#' GSEA v1 and is avoided by GSEA v2.  All results from FES or GSEA should be
#' evaluated by looking at the plots and changes in gene expression to evaluate
#' the feasibility of detecting the shift in the gene set distribution using
#' experimental tools.  Typically if the changes are smaller than a 0.5 on the
#' log2FC scale, they will be challenging to detect by methods such as Q-PCR.
#'
#'
#'
#' @param sn.table a table of ordered, valued gene rankings
#' @param geneset a set of genes to determine the FES score
#' @param hits.only Emit all genes(False) or only those which are in the geneset(True)
#'
#' @return A data.frame of FES scores for the selected genes.
#' @export
#'
#' @examples
fes.cmp <- function(sn.table, sn.table2) {
  order <- match(rownames(sn.table), rownames(sn.table2))
  hits <- rep(0, nrow(sn.table))
  hits2 <- rep(0, nrow(sn.table))
  lower <- rep(1, nrow(sn.table))
  higher <- rep(1, nrow(sn.table))
  for (index in 1:length(order)) {
    hits[index] <- sum(order[1:index] <= index)
    hits2[index] <- sum(order[index:nrow(sn.table)] <= nrow(sn.table) - index)
    lower[index] <- phyper(hits[index], nrow(sn.table), nrow(sn.table), index, lower.tail = F)
    higher[index] <- phyper(hits2[index], nrow(sn.table), nrow(sn.table), nrow(sn.table) - index, lower.tail = T)
  }
  data.frame(sn.table, p.lower = lower, p.higher = higher, hits = hits)
}
