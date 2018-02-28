# gseasier

## Overview
This R package contains everything you need to run [GSEA](https://software.broadinstitute.org/gsea/index.jsp) from inside of R.  GSEA stands for Gene Set Enrichment Analysis[1](http://www.pnas.org/content/102/43/15545) [2](https://www.nature.com/articles/ng1180).  These functions allow you to run GSEA from data you already have in R say from differential gene analysis.

## Installation
```
# The easiest method to install this package it from GitHub:
# install.packages("devtools")
devtools::install_github("jhart99/gseasier")
```
In addition you will need the GSEA java binary and signatures.  These must be obtained from The Broad Institute.

## Usage
The basic usage of gseasier is outlined below.  This details loading of genesets, data, calculation of a ranking, in this case using signal to noise, execution of GSEA and finally reading the results back into R.  Again the genesets are from Msigdb, the GSEA binaries are from The Broad and the example data is available from The Broad.
```
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
set.gsea("~/gsea/gsea-3.0.jar")
c1.gender  <- GSEA.preranked.java.exec(input.rnk = "gender.rnk",
                                     gs.db = "gsea/c1.all.v6.1.symbols.gmt",
                                     output.directory = "gsea.out",
                                     gs.size.threshold.max = 1500,
                                     doc.string = "c1.all.2.0")
# read in the data from the files
gender.gsea <- cbind(read.gsea.results(c1.gender),
                     cat = "c1.all.2.0", data = "gender")
                     
```

# Getting help

If there are problems with these routines, feel free to open an issue and I will try to solve it.  This is very much a work in progress and it is possible that functions might be added or removed or that function signatures may change in the future up to v1.0.
