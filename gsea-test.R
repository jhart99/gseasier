# Copyright ---------------------------------------------------------------
# 2014 The Scripps Research Institute Author: Jonathan Ross Hart

# Author ------------------------------------------------------------------
# Jonathan Ross Hart(jonathan@jonathanrosshart.com)

# Description -------------------------------------------------------------
# Fisher enrichment

# Input -------------------------------------------------------------------
# Gene expression data and genesets

# Methods -------------------------------------------------

# Outputs -------------------------------------------------


# Library imports ---------------------------------------------------------
library(gseasier)

# Load the MsigDB signatures
sig.gmts <- list.files("gsea", "all.*symbols")
sig.list <- lapply(sig.gmts, loadSig, path = "gsea")
sigs <- do.call(rbind, sig.list)
rm(sig.list, sig.gmts)

# Read gender data.  This is example data from Broad that is used for
# demonstration of the GSEA program.  You can get a copy here.
# http://software.broadinstitute.org/gsea/datasets.jsp
gender <- read.gct("datasets//Gender.gct")
gender.cls <- read.cls("datasets//Gender.cls")

# score it with signal to noise scores
gender.sn <- s2n.rank(gender, gender.cls)

# write the data to a file for running GSEA
write.sn.rnk(gender.sn, "gender.rnk")

# run GSEA against a signature database
set.gsea("~/gsea/gsea2-2.0.13.jar")
c1.gender  <- GSEA.preranked.java.exec(input.rnk = "gender.rnk",
                                     gs.db = "gsea/c1.all.v6.1.symbols.gmt",
                                     output.directory = "gsea.out",
                                     gs.size.threshold.max = 1500,
                                     doc.string = "c1.all.2.0")
# read in the data from the files
gender.gsea <- cbind(read.gsea.results(c1.gender),
                     cat = "c1.all.2.0", data = "gender")
set.gsea("~/gsea/gsea2-2.2.2.jar")
c1.gender  <- GSEA.preranked.java.exec(input.rnk = "gender.rnk",
                                       gs.db = "gsea/c1.all.v6.1.symbols.gmt",
                                       output.directory = "gsea.out",
                                       gs.size.threshold.max = 1500,
                                       doc.string = "c1.all.2.2")
# read in the data from the files
gender.gsea <- cbind(read.gsea.results(c1.gender),
                     cat = "c1.all.2.2", data = "gender")
set.gsea("~/gsea/gsea-3.0.jar")
c1.gender  <- GSEA.preranked.java.exec(input.rnk = "gender.rnk",
                                       gs.db = "gsea/c1.all.v6.1.symbols.gmt",
                                       output.directory = "gsea.out",
                                       gs.size.threshold.max = 1500,
                                       doc.string = "c1.all.3.0")
# read in the data from the files
gender.gsea <- cbind(read.gsea.results(c1.gender),
                     cat = "c1.all.3.0", data = "gender")

# make a plot of the top result.

# example methods for ES.plot  All of the following will produce the same plot
ES.plot(gender, gender.cls, NULL, sigs$gene[sigs$sig == "chr4q22"])
ES.plot(gender, gender.cls, geneset = sigs$gene[sigs$sig == "chr9p22"])
ES.plot(gender, gender.cls, geneset = "chr9p22", msigdb = sigs)

ES.plot(sn.table = gender.sn, geneset = sigs$gene[sigs$sig == "chr9p22"])
ES.plot(sn.table = gender.sn, geneset = "chr9p22", msigdb = sigs)

c1.gender  <- GSEA.preranked.java.exec(input.rnk = "gender.rnk",
                                       gs.db = "gsea/c2.cp.v6.1.symbols.gmt",
                                       output.directory = "gsea.out",
                                       gs.size.threshold.max = 1500,
                                       doc.string = "c1.all.3.0")
# read in the data from the files
gender.gsea <- cbind(read.gsea.results(c1.gender),
                     cat = "c2.cp", data = "gender")

# lets see if the top hits significantly overlap with each other
test.sigs <- sigs[sigs$sig %in% gender.gsea$sig[gender.gsea$q < 0.01], ]
# you need to know the set of all genes present in the signatures and your test
# set.  We need to know this because in the comparison stage the matrixes need
# the same dimensions and order.
genes <- make.names(unique(test.sigs$gene))
# convert the gene signature table into a sparse matrix
test.matrix <- convert.sigs.to.matrix(test.sigs, genes)
test.compare <- compare.sigs(test.matrix, test.matrix)

library(gplots)
library(reshape2)
test.comp.matrix <- acast(test.compare, sig ~ test.sig,
                          value.var = "q.value")
test.comp.matrix[test.comp.matrix == 0] <- .Machine$double.xmin
test.comp.matrix <- -log2(test.comp.matrix)
heatmap.2(test.comp.matrix, trace = "none", scale = "none",
          margins = c(12,12))

# Leading edge
sigs.to.test <- as.character(gender.gsea$sig[gender.gsea$q < 0.01])

leading.edge.heatmap(gender.sn, sigs[sigs$sig %in% sigs.to.test, ])

leading.edge(gender.sn, sigs$gene[sigs$sig == "REACTOME_METABOLISM_OF_RNA"])
leading.edge(gender.sn, sigs$gene[sigs$sig == "REACTOME_CELL_CYCLE"])


# FES FES is Fast Enrichment Scoring.  It uses an iterative Fisher Enrichment
# score over the preranked vector to look for enrichment of a signature's
# members near either end of the ranking.  FES is not the same algorithm as
# GSEA, but can be a fast way to look at all of the signatures prior to doing
# GSEA.  You might also filter the signatures you look at by GSEA using this
# method.

gender.fes.table <- FES.table(gender.sn, sigs[sigs$class == "c1.all",])

