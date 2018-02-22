# Copyright ---------------------------------------------------------------
# 2016 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------

# Function to write a new gmt file

# Input -------------------------------------------------------------------

# a table of signatures from msgidb

# Methods -------------------------------------------------

# write.gmt

# Outputs -------------------------------------------------

# Returns a vector of factors deciding if the gene is up, down or unchanged.
# The values are multiple test corrected.

#' Title
#'
#' @param in.sigs
#' @param out.file
#'
#' @return
#' @export
#'
#' @examples
write.gmt <- function(in.sigs, out.file) {
  library(plyr)
  colnames(in.sigs) <- c("sig", "sig.class", "def", "gene")
  file.handle <- file(out.file, "wt")
  plyr::dlply(in.sigs, .(sig.class, sig), .fun = function(df) {
    write(c(paste0(df$sig.class[1], ".", df$sig[1]),
            as.character(df$def[1]),
            as.character(unlist(df$gene))),
          file = file.handle, append = T, sep = "\t", ncolumns = 5000)
  })
  close(file.handle)
}
