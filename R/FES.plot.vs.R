# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A gene set overlap function


#' Title
#'
#' @param geneset
#' @param sn.table.list
#'
#' @return
#' @export
#'
#' @examples
FES.plot.vs <- function(sn.table.list, geneset) {
  plot.data <- do.call(
    rbind, lapply(seq_along(sn.table.list),
                  function(d) {
                    out <- FES(sn.table.list[[d]], geneset, hits.only = F)
                    out$sn.table <- d
                    out
                  }))
  ES.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = rank)) +
    ggplot2::geom_line(color = "blue", aes(y = -log10(p.lower))) +
    ggplot2::geom_line(color = "red", aes(y = -log10(p.higher))) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~sn.table) +
    ggplot2::ylab(expression(log[10]*p))
  ES.plot
}

