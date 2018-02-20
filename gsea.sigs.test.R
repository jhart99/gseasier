# Copyright ---------------------------------------------------------------
# 2018 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# A demonstration of the gene signature overlap routines. 

# Input -------------------------------------------------------------------
# msigdb formatted gmt files or yuor own gmt formatted gene sets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------

# a table of comparisons between two sets of gene sets with odds ratios and
# p.values

# Library imports ---------------------------------------------------------

# acquire the msigdb files from Broad directly at
# http://software.broadinstitute.org/gsea/downloads.jsp

# this part will work with either set of genes.  Meaning you can use entrez ids
# or gene names as you prefer.

source("R/msigdb.R")
sig.gmts <- list.files("gsea", "all.*entrez")
sig.list <- lapply(sig.gmts, loadSig, path = "gsea")
sigs <- do.call(rbind, sig.list)
rm(sig.list, sig.gmts)

source("R/gsea.sigs.R")
# I'm selecting a list of signatures that is just the set of all the hallmarks.
test.sigs <- sigs[sigs$class == "h.all", ]
# you need to know the set of all genes present in the signatures and your test
# set.  We need to know this because in the comparison stage the matrixes need
# the same dimensions and order.
genes <- make.names(unique(c(sigs$gene, test.sigs$gene)))
# convert the gene signature table into a sparse matrix
sig.matrix <- convert.sigs.to.matrix(sigs, genes)
test.matrix <- convert.sigs.to.matrix(test.sigs, genes)

# make the comparison and get the odds ratios and p.values
test.compare <- compare.sigs(sig.matrix, test.matrix)

head(test.compare[order(test.compare$q.value), ], 75)

toot <- test.matrix %*% t(test.matrix)
