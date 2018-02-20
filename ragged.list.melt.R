# Copyright ---------------------------------------------------------------
# 2013-2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Melts a ragged list into a long matrix

# Input -------------------------------------------------------------------
# Ragged list object

# Methods -------------------------------------------------
# ragged.list.melt

# Outputs -------------------------------------------------
# a matrix/data.frame which is long and with a fixed number of columns



require(foreach)
require(doMC)
registerDoMC(cores = 8)
ragged.list.melt <- function(x) {
  do.call(rbind,
          lapply(1:length(x), function(i) {
            if (length(x[[i]]) > 0) {
              cbind(i, x[[i]])
            }
          }))
}