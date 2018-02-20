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

# Load the MsigDB signatures
source("R/msigdb.R")
sig.gmts <- list.files("gsea", "all.*symbols")
sig.list <- lapply(sig.gmts, loadSig, path = "gsea")
sigs <- do.call(rbind, sig.list)
rm(sig.list, sig.gmts)

# Read gender data.  This is example data from Broad that is used for
# demonstration of the GSEA program.  You can get a copy here.
# http://software.broadinstitute.org/gsea/datasets.jsp
source("R/read.gct.R")
source("R/read.cls.R")
gender <- read.gct("datasets//Gender.gct")
gender.cls <- read.cls("datasets//Gender.cls")

# score it with signal to noise scores
source("R/s2n.rank.R")
source("R/write.rnk.R")
source("R/gsea.java.R")
gender.sn <- s2n.rank(gender, gender.cls)

# write the data to a file for running GSEA
write.sn.rnk(gender.sn, "gender.rnk")

# run GSEA against a signature database
c1.gender  <- GSEA.preranked.java.exec(input.rnk = "gender.rnk",
                                     gs.db = "gsea/c1.all.v6.1.symbols.gmt",
                                     output.directory = "gsea.out",
                                     gs.size.threshold.max = 1500,
                                     doc.string = "c1.all")

# read in the data from the files
gender.gsea <- cbind(read.gsea.results(c1.gender), 
                     cat = "c1.all", data = "gender")

# make a plot of the top result.
source('R/ES.plot.R')
ES.plot(gender, gender.cls, NULL, sigs$gene[sigs$sig == "chr9p22"])
