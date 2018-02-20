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

write.gmt <- function(in.sigs, out.file) {
  library(plyr)
  colnames(in.sigs) <- c("sig", "sig.class", "def", "gene")
  file.handle <- file(out.file, "wt")
  dlply(in.sigs, .(sig.class, sig), .fun=function(df) {
    write(c(paste0(df$sig.class[1], ".", df$sig[1]), as.character(df$def[1]), as.character(unlist(df$gene))), file=file.handle, append=T, sep="\t", ncolumns=5000)
    #paste(paste0(class,df$sig[1]), df$def[1], df$gene,sep="\t")
  })
  close(file.handle)
}