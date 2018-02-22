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
ES.plot(gender, gender.cls, NULL, sigs$gene[sigs$sig == "chr9p22"])
ES.plot(gender, gender.cls, geneset = sigs$gene[sigs$sig == "chr9p22"])
ES.plot(gender, gender.cls, geneset = "chr9p22", msigdb = sigs)

ES.plot(sn.table = gender.sn, geneset = sigs$gene[sigs$sig == "chr9p22"])
ES.plot(sn.table = gender.sn, geneset = "chr9p22", msigdb = sigs)
