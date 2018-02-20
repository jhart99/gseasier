write.gct <- function(in.DGE, out.file) {
  file.handle <- file(out.file, "wt")
  writeLines("#1.2", file.handle)
  cat(c(dim(log(in.DGE$pseudo.counts)), "\n"), file=file.handle)
  gmt.table <- data.frame(NAME=rownames(in.DGE$pseudo.counts), DESC="(na)", in.DGE$pseudo.counts)
  
  write.table(gmt.table, file=file.handle, sep="\t", row.names=F, quote=F)
  close(file.handle)
}