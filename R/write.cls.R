write.cls <- function(in.DGE, out.file) {
  file.handle <- file(out.file, "wt")
  cat(c(length(RNA.DGE$samples$group), length(levels(RNA.DGE$samples$group)), "1\n"), file=file.handle)
  cat(c("#", levels(RNA.DGE$samples$group), "\n"), file=file.handle)
  cat(as.character(RNA.DGE$samples$group), file=file.handle)
  close(file.handle)
}