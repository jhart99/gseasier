# Copyright ---------------------------------------------------------------
# 2013-2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Write a rank file from expression data for preranked GSEA

# Input -------------------------------------------------------------------
# expression.data 

# Methods -------------------------------------------------
# write.rnk

# Outputs -------------------------------------------------
# rank file

write.rnk <- function(expression.data, out.file, type) {
  file.handle <- file(out.file, "wt")
  writeLines("#", file.handle)
  #cat(c(dim(log(in.DGE$pseudo.counts)), "\n"), file=file.handle)
  if(type == "prot") rnk.table <- data.frame(gene=expression.data$gene.x, expr=expression.data$median)
  if(type == "rna") rnk.table <- data.frame(gene=expression.data$gene, expr=expression.data$rna.sn)
  rnk.table <- rnk.table[order(rnk.table$expr), ]
  rnk.table <- rnk.table[!is.na(rnk.table$expr), ]
  write.table(rnk.table, file=file.handle, sep="\t", row.names=F, quote=F)
  close(file.handle)
}
write.sn.rnk <- function(expression.data, out.file) {
  file.handle <- file(out.file, "wt")
  writeLines("#", file.handle)
  rnk.table <- data.frame(gene=rownames(expression.data), expr=expression.data$sn)
  rnk.table <- rnk.table[order(rnk.table$expr), ]
  rnk.table <- rnk.table[!is.na(rnk.table$expr), ]
  write.table(rnk.table, file=file.handle, sep="\t", row.names=F, quote=F)
  close(file.handle)
}