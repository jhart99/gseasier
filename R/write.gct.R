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
write.gct <- function(in.DGE, out.file) {
  file.handle <- file(out.file, "wt")
  writeLines("#1.2", file.handle)
  cat(c(dim(log(in.DGE$pseudo.counts)), "\n"), file = file.handle)
  gmt.table <- data.frame(NAME = rownames(in.DGE$pseudo.counts),
                          DESC = "(na)", in.DGE$pseudo.counts)

  write.table(gmt.table, file = file.handle, sep = "\t",
              row.names = F, quote = F)
  close(file.handle)
}
