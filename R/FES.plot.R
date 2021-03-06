# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A demonstration of the gene signature overlap routines.


#' FES plot
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
#' @param gmt a gmt formatted gene expression data.frame
#' @param cls a vector describing the classes the gmt datasets belong to
#' @param comparison which comparison of classes is desired defaults to the
#'   first two levels from cls
#' @param geneset the set of genes to test and plot
#' @param sn.table the ordered, valued ranking of genes
#' @param msigdb a database of signatures
#'
#' @return A ggplot2 plot of the FES score in a style similar to the GSEA style
#'   plots
#' @export
#'
#' @examples
FES.plot <- function(gmt=NULL, cls=NULL, comparison=NULL,
                     geneset=NULL, sn.table=NULL, msigdb=NULL) {
  if (is.null(gmt)) {
    if (is.null(sn.table)) {
      stop("Either gmt and cls or sn.table must be used to generate a plot")
    }
  } else {
    if (length(comparison) == 2 | is.null(comparison)) {
      sn.table <- s2n.rank(gmt, cls, comparison)
    } else {
      stop("Invalid comparison")
    }
  }

  if (length(geneset) == 1) {
    if (is.null(msigdb)) {
      stop("msigdb missing")
    }
    geneset <- msigdb$gene[msigdb$sig == geneset]
  }

  plot.data <- FES(sn.table, geneset, hits.only = F)

  plot.data$percentile <- cut(plot.data$sn,
                              quantile(plot.data$sn, probs = seq(0, 1, 1 / 9)))

  grid::grid.newpage()
  theme.blank.x <- ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(0, 1, 0, 0), "lines"),
    panel.border = ggplot2::element_rect(size = 0),
    panel.grid = ggplot2::element_blank())
  theme.blank.y <- theme.blank.x +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  sn.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = rank, y = sn)) +
    ggplot2::geom_point() + ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 1, 0, 0), "lines"),
                   panel.border = ggplot2::element_rect(size = 0),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::ylab("Signal to Noise")
  ES.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = rank)) +
    ggplot2::geom_line(color = "blue", aes(y = -log10(p.lower))) +
    ggplot2::geom_line(color = "red", aes(y = -log10(p.higher))) +
    ggplot2::theme_bw() +
    theme.blank.x +
    ggplot2::ylab(expression(log[10]*p))
  ticks <- ggplot2::ggplot(plot.data, ggplot2::aes(x = rank,
                                                   xend = rank,
                                                   y = ifelse(hit, 0.1, 0),
                                                   yend = 0)) +
    ggplot2::geom_segment() + ggplot2::theme_bw() +
    theme.blank.y
  gradient <- ggplot2::ggplot(plot.data, ggplot2::aes(x = rank)) +
    ggplot2::theme_bw() +
    ggplot2::geom_tile(ggplot2::aes(y = 0, fill = as.numeric(percentile))) +
    ggplot2::scale_fill_gradient2(midpoint = 5, low = "blue", high = "red") +
    theme.blank.y +
    ggplot2::theme(legend.position = "none")

  # There are issues getting the 4 plots to line up nicely.  This is one way to
  # make it work.  It may break and there might be a better way to do this, but
  # this is what I have.  What is going on is we make grobs of our plots then
  # find the widths for the plot areas.  Then set all of the plot areas to have
  # the same widths.  When these are plotted the 4 graphs will line up.

  gsn <- ggplot2::ggplotGrob(sn.plot)
  gES <- ggplot2::ggplotGrob(ES.plot)
  gticks <- ggplot2::ggplotGrob(ticks)
  ggradient <- ggplot2::ggplotGrob(gradient)
  groblist <- list(gES, gticks, ggradient, gsn)
  widths <- list()
  for (i in 1:length(groblist)) {
    widths[[i]] <- groblist[[i]]$widths[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  for (i in 1:length(groblist)) {
    groblist[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  grid::grid.newpage()
  # why is this manually unrolled?  Well grid looks at groblist and decides it
  # is a gList so there is no way for me to pass in the other parameters for
  # grid.arrange because it will error and say these are grobs in the gList.
  gridExtra::grid.arrange(groblist[[1]], groblist[[2]], groblist[[3]],
                          groblist[[4]], ncol = 1, nrow = 4,
                          heights = c(5, 1, 1, 3))

}
