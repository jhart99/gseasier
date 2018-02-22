# Copyright ---------------------------------------------------------------
# 2013-2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Melts a ragged list into a long matrix

# Input -------------------------------------------------------------------
# Ragged list object of constant type(numeric, character, etc)

# Methods -------------------------------------------------
# ragged.list.melt

# Outputs -------------------------------------------------
# a

#' Title
#'
#' @param x A list of uneven vectors of constant type
#'
#' @return
#' @export
#'
#' @examples
ragged.list.melt <- function(x) {
  do.call(rbind,
          lapply(seq_along(x), function(i) {
            if (length(x[[i]]) > 0) {
              cbind(i, x[[i]])
            }
          }))
}
