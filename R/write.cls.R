# Copyright ---------------------------------------------------------------
# 2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# signal to noise scoring function

# Input -------------------------------------------------------------------
# Gene expression data and genesets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------


# Library imports ---------------------------------------------------------

# Determine signal to noise scoring for each line in the data table.  dataA is
# the baseline(negative if higher) and dataB is the test(positive if higher)

#' Title
#'
#' @param in.DGE
#' @param out.file
#'
#' @return
#' @export
#'
#' @examples
write.cls <- function(in.DGE, out.file) {
  file.handle <- file(out.file, "wt")
  cat(c(length(in.DGE$samples$group),
        length(levels(in.DGE$samples$group)), "1\n"), file = file.handle)
  cat(c("#", levels(in.DGE$samples$group), "\n"), file = file.handle)
  cat(as.character(in.DGE$samples$group), file = file.handle)
  close(file.handle)
}
