# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Enrichment set plot similar to the one provided by the GSEA program

#' Enrichment set plot
#'
#' @param gmt A gmt formatted gene expression data set
#' @param cls A cls formatted classifier
#' @param comparison Which levels of cls should be compared defaults to levels 1 and 2
#' @param geneset Either a vector of gene identifiers or a geneset name in combination with msigdb parameter
#' @param sn.table (optional) A signal to noise named vector can be used in place of gmt, cls and comparison
#' @param msigdb (optional) A database of signatures
#'
#' @return A ggplot2 plot of the GSEA plot
#' @export
#'
#' @examples
#' ES.plot(gender, gender.cls, NULL, sigs$gene[sigs$sig == "chr4q22"])
#' ES.plot(gender, gender.cls, geneset = sigs$gene[sigs$sig == "chr9p22"])
#' ES.plot(gender, gender.cls, geneset = "chr9p22", msigdb = sigs)
#'
ES.plot <- function(gmt=NULL, cls=NULL, comparison=NULL,
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

  plot.data <- ES(sn.table, geneset, hits.only = F)

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
          panel.grid = ggplot2::element_blank())
  ES.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = rank, y = ES)) +
    ggplot2::geom_line(color = "blue") + ggplot2::theme_bw() +
    theme.blank.x
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

    # hacky things to modify the facets to different heights

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
