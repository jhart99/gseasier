# Copyright ---------------------------------------------------------------
# 2013 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Function for reading msigdb gmt files

# Input -------------------------------------------------------------------
# msigdb formatted gmt file

# Methods -------------------------------------------------

# Outputs -------------------------------------------------
# returns a data frame of sigs with genes



# Library imports ---------------------------------------------------------

loadSig <- function(filename, path=".") {
  source("ragged.list.melt.R")
  
  sig.file <- file(paste0(path, "/", filename), "rt")
  sig.text <- readLines(sig.file)
  close(sig.file)
  
  split.lines <- strsplit(sig.text, "\t")
  sig <- unlist(lapply(split.lines, "[", 1))
  definitions <- unlist(lapply(split.lines, "[", 2))
  genes <- lapply(split.lines, "[", c(-1,-2))
  
  genes <- data.frame(ragged.list.melt(genes), stringsAsFactors=F)
  genes <- data.frame(lapply(genes, type.convert, as.is=T))
  genes$sig <- sig[genes$i]
  genes$def <- definitions[genes$i]
  filename <- paste(unlist(strsplit(filename, "\\."))[c(1,2)], collapse=".")
  data.frame(sig=genes$sig, class=filename, def=genes$def, gene=genes$cur.line)
}

plotSigOverlap <- function(sigs) {
  library(foreach)
  library(gplots)
  sigList <- as.character(unique(sigs$sig))
  graph.data <- foreach::foreach(sig1=iter(sigList), .combine=rbind) %:%
    foreach::foreach(sig2=iter(sigList), .combine=c) %do% {
      total <- sum(sigs$sig %in% c(sig1, sig2))
      shared <- length(unique(sigs$gene[sigs$sig %in% c(sig1, sig2)]))
      mine <- sum(sigs$sig %in% c(sig1))
      1-(total-shared)/mine
    }
  rownames(graph.data) <- sigList
  colnames(graph.data) <- sigList
  heatmap.2(graph.data, trace="none")
  
}