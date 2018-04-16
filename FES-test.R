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
library(edgeR)
library(reshape)
library(plyr)
library(foreach)
library(doMC)
registerDoMC(cores=8)
library(ggplot2)
library(hexbin)
library(RColorBrewer)
library(boot)
library(doBy)
library(Rcpp)
library(gplots)

source("GSEA.java.R")
source("msigdb.R")
source("ragged.list.melt.R")
source("FES.R")
source("read.gct.R")
source("read.cls.R")
source("read.rnk.R")
source("NES.R")
source("ES.R")
source("tcga.add.maf.R")
source("read.tcga.cls.R")
# Load the MsigDB signatures
sigs <- loadSig("c2.cp.v4.0.symbols.gmt")
sigs <- rbind(sigs, loadSig("c2.cgp.v4.0.symbols.gmt"))
sigs <- rbind(sigs, loadSig("c1.all.v4.0.symbols.gmt"))
sigs <- rbind(sigs, loadSig("c3.all.v4.0.symbols.gmt"))
sigs <- rbind(sigs, loadSig("c5.all.v4.0.symbols.gmt"))
sigs <- rbind(sigs, loadSig("c6.all.v4.0.symbols.gmt"))

# Read gender data
gender <- read.gct("datasets//Gender.gct")
gender.cls <- read.cls("datasets//Gender.cls")

gender.sn <- s2n.rank(gender, gender.cls)

gender.table <- FES.table(gender.sn, sigs[sigs$class="c1.all.v4.0.symbols.gmt",])
head(gender.table[order(gender.table$p.hi),])
head(gender.table[order(gender.table$p.lo),])
FES.plot2(gender.sn, sigs$gene[sigs$sig=="chryq11"])
FES.plot(gender.sn, sigs$gene[sigs$sig=="chryp11"])
FES.plot(gender.sn, sigs$gene[sigs$sig=="SU_TESTIS"])

FES.plot(gender.sn, sigs$gene[sigs$sig=="DISTECHE_ESCAPED_FROM_X_INACTIVATION"])
FES.plot(gender.sn, sigs$gene[sigs$sig=="RUNNE_GENDER_EFFECT_UP"])

# Read Leukemia data
leukemia <- read.gct("datasets//Leukemia.gct")
leukemia.cls <- read.cls("datasets//Leukemia.cls")

leukemia.sn <- s2n.rank(leukemia, leukemia.cls)

ggplot(leukemia.sn, aes(x=rank, y=sn)) + geom_point()

leukemia.table <- FES.table(leukemia.sn, sigs[sigs$class=="c1.all.v4.0.symbols.gmt",])
head(leukemia.table[order(leukemia.table$p.hi),])
head(leukemia.table[order(leukemia.table$p.lo),])
FES.plot(leukemia.sn, sigs$gene[sigs$sig=="chr6q21"])
FES.plot(leukemia.sn, sigs$gene[sigs$sig=="chr5q31"])
FES.plot(leukemia.sn, sigs$gene[sigs$sig=="chr6p22"])
FES.plot(leukemia.sn, sigs$gene[sigs$sig=="chr19p13"])
FES.plot(leukemia.sn, sigs$gene[sigs$sig=="chr19q13"])
FES.plot(leukemia.sn, sigs$gene[sigs$sig=="chr16q13"])
# Read p53 data
p53 <- read.gct("datasets//P53.gct")
p53.cls <- read.cls("datasets//P53.cls")

p53.sn <- s2n.rank(p53, p53.cls)

ggplot(p53.sn, aes(x=rank, y=sn)) + geom_point()

p53.table <- FES.table(p53.sn, sigs[sigs$class=="c2.cp.v4.0.symbols.gmt",])
head(p53.table[order(p53.table$p.hi),],20)
head(p53.table[order(p53.table$p.lo),], 20)
FES.plot(p53.sn, sigs$gene[sigs$sig=="P53_DN.V1_UP"])
FES.plot(p53.sn, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])

png("test.png", 1500, 1500, res=300)
FES.plot(p53.sn, sigs$gene[sigs$sig=="P53_DN.V1_DN"])
dev.off()
FES.plot(p53.sn, sigs$gene[sigs$sig=="REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS"])
FES.plot(p53.sn, sigs$gene[sigs$sig=="PID_P53DOWNSTREAMPATHWAY"])

# Read Lung from Boston data
lung.boston <- read.gct("datasets//Lung_Boston.gct")
lung.boston.cls <- read.cls("datasets//Lung_Boston.cls")

lung.boston.sn <- s2n.rank(lung.boston, lung.boston.cls)

ggplot(lung.boston.sn, aes(x=rank, y=sn)) + geom_point()

lung.boston.table <- FES.table(lung.boston.sn, sigs)
head(lung.boston.table[order(lung.boston.table$p.hi),])
head(lung.boston.table[order(lung.boston.table$p.lo),])
FES.plot(lung.boston.sn, sigs$gene[sigs$sig=="BIOCARTA_VEGF_PATHWAY"])

# Read Lung from Michigan data
lung.michigan <- read.gct("datasets//Lung_Michigan.gct")
lung.michigan.cls <- read.cls("datasets//Lung_Michigan.cls")

lung.michigan.sn <- s2n.rank(lung.michigan, lung.michigan.cls)

ggplot(lung.michigan.sn, aes(x=rank, y=sn)) + geom_point()

lung.michigan.table <- FES.table(lung.michigan.sn, sigs)
head(lung.michigan.table[order(lung.michigan.table$p.hi),])
head(lung.michigan.table[order(lung.michigan.table$p.lo),])
FES.plot(lung.michigan.sn, sigs$gene[sigs$sig=="KEGG_GLYCOLYSIS_GLUCONEOGENESIS"])

# P4936 data
p4936 <- read.gct("datasets//p4936.gct")
p4936.cls <- read.cls("datasets//p4936.cls")
p4936.sn.kj9 <- s2n.rank(p4936, p4936.cls, c("KJ9", "None"))

ggplot(p4936.sn.kj9, aes(x=rank, y=sn)) + geom_point()

p4936.kj9.table <- FES.table(p4936.sn.kj9, sigs)
p4936.kj9.table <- FES.table(p4936.sn.kj9, sigs[grep("MYC_",sigs$sig),])
head(p4936.kj9.table[order(p4936.kj9.table$p.hi),], 100)
head(p4936.kj9.table[order(p4936.kj9.table$p.lo),], 100)
png("test.png", 1500, 1500, res=300)
FES.plot(p4936.sn.kj9, sigs$gene[sigs$sig=="DANG_MYC_TARGETS_UP"])
dev.off()
FES.plot(p4936.sn.kj9, sigs$gene[sigs$sig=="SCHUHMACHER_MYC_TARGETS_UP"])
FES.plot(p4936.sn.kj9, sigs$gene[sigs$sig=="WEI_MYCN_TARGETS_WITH_E_BOX"])
FES.plot(p4936.sn.kj9, sigs$gene[sigs$sig=="ODONNELL_TARGETS_OF_MYC_AND_TFRC_DN"])
FES.plot(p4936.sn.kj9, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])

p4936.sn.dox <- s2n.rank(p4936, p4936.cls, c("DOX", "None"))
ggplot(p4936.sn.dox, aes(x=rank, y=sn)) + geom_point()
p4936.dox.table <- FES.table(p4936.sn.dox, sigs[grep("MYC_",sigs$sig),])
head(p4936.dox.table[order(p4936.dox.table$p.hi),], 100)
head(p4936.dox.table[order(p4936.dox.table$p.lo),], 100)
png("test.png", 1500, 1500, res=300)
FES.plot(p4936.sn.dox, sigs$gene[sigs$sig=="DANG_MYC_TARGETS_UP"])
dev.off()
FES.plot(p4936.sn.dox, sigs$gene[sigs$sig=="SCHUHMACHER_MYC_TARGETS_UP"])
FES.plot(p4936.sn.dox, sigs$gene[sigs$sig=="WEI_MYCN_TARGETS_WITH_E_BOX"])
FES.plot(p4936.sn.dox, sigs$gene[sigs$sig=="ODONNELL_TARGETS_OF_MYC_AND_TFRC_DN"])
FES.plot(p4936.sn.dox, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])


mtor0.silac.sn <- read.rnk("datasets//silac.combined.rnk")
ggplot(mtor0.silac.sn, aes(x=rank, y=sn)) + geom_point()
mtor0.silac.table <- FES.table(mtor0.silac.sn, sigs)
head(mtor0.silac.table[order(mtor0.silac.table$p.hi),], 20)
head(mtor0.silac.table[order(mtor0.silac.table$p.lo),], 20)
sum(mtor0.silac.table$q.hi<0.01)
sum(mtor0.silac.table$q.lo<0.01)
FES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN"])

png("charafe.silac.png", 4,5,units="in",res=300)
ES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN"])
dev.off()

FES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="NUCLEUS"])
FES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="CYTOPLASMIC_PART"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
FES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
NES(mtor0.rna.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
NES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
NES(mtor0.rna.sn, sigs$gene[sigs$sig=="REACTOME_CHROMOSOME_MAINTENANCE"])
NES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="REACTOME_CHROMOSOME_MAINTENANCE"])
mtor0.silac.table.NES <- NES.table(mtor0.silac.sn, sigs)
head(mtor0.silac.table.NES[order(mtor0.silac.table.NES$p),], 20)
sum(mtor0.silac.table.NES$q<0.01)

mtor0.rna <- read.gct("datasets//hvw.gct")
# this block median scales the data after log2 transforming it.  
# I am performing this only to show that given the same scales that there is no difference in the results.
mtor0.rna.scaled <- mtor0.rna
mtor0.rna.scaled[, c(-1, -2)] <- log2(mtor0.rna.scaled[,3:8])
shift <- apply(mtor0.rna.scaled[, c(-1, -2)], 1, median)
mtor0.rna.scaled <- foreach(i=1:nrow(mtor0.rna.scaled), .combine=rbind) %do% {
  mtor0.rna.scaled[i, c(-1, -2)] - shift[i]
}
mtor0.rna.scaled <- data.frame(mtor0.rna[, c(1, 2)], mtor0.rna.scaled)

mtor0.rna.cls <- read.cls("datasets//hvw.cls")
mtor0.rna.sn.scaled <- s2n.rank(mtor0.rna.scaled, mtor0.rna.cls)
mtor0.rna.table.scaled <- FES.table(mtor0.rna.sn.scaled, sigs)
head(mtor0.rna.table.scaled[order(mtor0.rna.table.scaled$p.hi),], 20)
head(mtor0.rna.table.scaled[order(mtor0.rna.table.scaled$p.lo),], 20)

foo.sn <- data.frame(sn=mtor0.rna.scaled[, 3], row.names=mtor0.rna.scaled[, 1])
foo.sn$rank <- 1:nrow(foo.sn)
foo.sn <- foo.sn[order(-foo.sn$sn), ]
foo.sn$rank <- 1:nrow(foo.sn)
foo <- FES.heatmap(tcga.breast, sigs[grep("SMID", sigs$sig), ])
row.cl <- hclust(dist(foo))
col.cl <- hclust(dist(t(foo)))
foo.col <- colnames(foo)
foo.cls <- as.factor(tcga.txt$PAM50[gsub("\\.", "-", colnames(foo)) %in% tcga.txt$Sample])
foo.colors <- brewer.pal(5, "Set1")
heatmap.2(foo[, c(1:10)], trace="none", col=colorRampPalette(c("green", "black", "red"))(n = 20), 
          breaks=-10:10, margins=c(5,20),
          Colv=FALSE, Rowv=FALSE, dendrogram="none")

ES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN"])


png("charafe.rna.png", 4,5,units="in",res=300)
ES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN"])
dev.off()

sum(mtor0.rna.table$q.hi<0.01)
sum(mtor0.rna.table$q.lo<0.01)
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])
foo <- FES.leading.edge(mtor0.rna.sn, sigs$gene[sigs$sig=="HUPER_BREAST_BASAL_VS_LUMINAL_UP"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="HUPER_BREAST_BASAL_VS_LUMINAL_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="HUPER_BREAST_BASAL_VS_LUMINAL_DN"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_DN"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="MYC_UP.V1_DN"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="ODONNELL_TARGETS_OF_MYC_AND_TFRC_UP"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="ODONNELL_TARGETS_OF_MYC_AND_TFRC_DN"])

# TCGA breast data
tcga.breast <- read.tcga("BRCA.exp.547.med.txt")
tcga.cls <- read.tcga.cls("BRCA.547.PAM50.SigClust.Subtypes.txt", "PAM50", tcga.breast)
tcga.maf <- read.maf("genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf")
tcga.txt <- read.tcga.txt("BRCA.547.PAM50.SigClust.Subtypes.txt")

# Basal will be positive(high) Luminal will be negative(low)
tcga.bvl.sn <- s2n.rank(tcga.breast, tcga.cls, c("Basal", "LumA"))
tcga.bvl.table <- FES.table(tcga.sn, sigs)
head(tcga.bvl.table[order(tcga.table$p.hi),], 20)
head(tcga.bvl.table[order(tcga.table$p.lo),], 20)

# Basal will be positive(high) normal will be negative(low)
tcga.bvn.sn <- s2n.rank(tcga.breast, tcga.cls, c("Basal", "Normal"))
tcga.bvn.table <- FES.table(tcga.bvn.sn, sigs)
head(tcga.bvn.table[order(tcga.bvn.table$p.hi),], 20)
head(tcga.bvn.table[order(tcga.bvn.table$p.lo),], 20)


# giant dataset
giant <- read.csv("giant.csv", quote="\"")
giant <- data.frame(NAME=giant$X, DESC="", giant[, -1])
giant.cls <- factor(c("h0n", "h0n", "h0n",
                      "h12a", "h12a", "h12a",
                      "h12n", "h12n", "h12n",
                      "h12r", "h12r", "h12r",
                      "h24a", "h24a", "h24a",
                      "h24n", "h24n", "h24n",
                      "h24r", "h24r", "h24r",
                      "h6a", "h6a", "h6a",
                      "h6n", "h6n", "h6n",
                      "h6r", "h6r", "h6r",
                      "w0n", "w0n", "w0n",
                      "w12a", "w12a", "w12a",
                      "w12n", "w12n", "w12n",
                      "w12r", "w12r", "w12r",
                      "w24a", "w24a", "w24a",
                      "w24n", "w24n", "w24n",
                      "w24r", "w24r", "w24r",
                      "w6a", "w6a", "w6a",
                      "w6n", "w6n", "w6n",
                      "w6r", "w6r", "w6r"
                      ))
giant.scaled <- giant

giant.scaled <- giant.scaled[apply(giant.scaled[,c(-1, -2)], 1, min) > 1, ] # We need to do this to eliminate a few which have a floating point 0 which kills the scaling

giant.scaled[, c(-1, -2)] <- log2(giant.scaled[,c(-1, -2)])
# is.na(giant.scaled) <- sapply(giant.scaled, is.infinite)
shift <- apply(giant.scaled[, c(-1, -2)], 1, median)
giant.scaled.names <- giant.scaled[, c(1, 2)]
giant.scaled <- foreach(i=1:nrow(giant.scaled), .combine=rbind) %dopar% {
  giant.scaled[i, c(-1, -2)] - shift[i]
}

giant.scaled <- data.frame(giant.scaled.names[, c(1, 2)], giant.scaled)
rownames(giant.scaled) <- giant.scaled$NAME

# giant.scaled.sn <- s2n.rank(giant.scaled, giant.cls, list(c("h0n", "h6n", "h12n", "h24n"), c("w0n", "w6n", "w12n", "w24n")))
# giant.scaled.table <- FES.table(giant.scaled.sn, sigs)
# heatmap.2(as.matrix(giant.scaled[giant.scaled$NAME %in% sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"], c(3:5, 27:29, 9:11, 18:20)]), 
#            trace="none", scale="none", col=colorRampPalette(c("blue", "white", "red")), labRow="", key=F, Colv=F)
# 
# heatmap.2(as.matrix(giant.scaled[giant.scaled$NAME %in% sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_DN"], c(3:5, 27:29, 9:11, 18:20)]), 
#           breaks=seq(-5,5,0.1), trace="none", scale="none", col=colorRampPalette(c("blue", "white", "red")), labRow="", key=F, Colv=F)
# 
# heatmap.2(as.matrix(giant.scaled[giant.scaled$NAME %in% sigs$gene[sigs$sig %in% c("SMID_BREAST_CANCER_BASAL_UP")], c(3:5, 9:11, 18:20, 27:29, 33:35, 39:41, 48:50, 57:60)]), 
#           breaks=seq(-3,3,0.1), trace="none", scale="none", col=colorRampPalette(c("blue", "white", "red")))
# 
# giant.scaled.early.v.late.sn <- s2n.rank(giant.scaled, giant.cls, list(c("h0n", "h6n", "w0n", "w6n"), c("h12n", "h24n", "w12n", "w24n")))
# giant.scaled.early.v.late.table <- FES.table(giant.scaled.early.v.late.sn, sigs)
# head(giant.scaled.early.v.late.table[order(giant.scaled.early.v.late.table$p.lo), ], 20)
# head(giant.scaled.early.v.late.table[order(giant.scaled.early.v.late.table$p.hi), ], 20)

giant.0.sn <- s2n.rank(giant.scaled, giant.cls, list("h0n", "w0n"))
giant.6.sn <- s2n.rank(giant.scaled, giant.cls, list("h6n", "w6n"))
giant.12.sn <- s2n.rank(giant.scaled, giant.cls, list("h12n", "w12n"))
giant.24.sn <- s2n.rank(giant.scaled, giant.cls, list("h24n", "w24n"))
giant.0.table <- FES.table(giant.0.sn, sigs)
giant.6.table <- FES.table(giant.6.sn, sigs)
giant.12.table <- FES.table(giant.12.sn, sigs)
giant.24.table <- FES.table(giant.24.sn, sigs)

giant.table <- data.frame(sig=giant.0.table$sig)
giant.table <- cbind(giant.table, h0.lo=log(giant.0.table$p.lo), h0.hi=-log(giant.0.table$p.hi))
giant.table <- cbind(giant.table, h6.lo=log(giant.6.table$p.lo), h6.hi=-log(giant.6.table$p.hi))
giant.table <- cbind(giant.table, h12.lo=log(giant.12.table$p.lo), h12.hi=-log(giant.12.table$p.hi))
giant.table <- cbind(giant.table, h24.lo=log(giant.24.table$p.lo), h24.hi=-log(giant.24.table$p.hi))
giant.table$h0 <- ifelse(-giant.table$h0.lo>giant.table$h0.hi, giant.table$h0.lo, giant.table$h0.hi)
giant.table$h6 <- ifelse(-giant.table$h6.lo>giant.table$h6.hi, giant.table$h6.lo, giant.table$h6.hi)
giant.table$h12 <- ifelse(-giant.table$h12.lo>giant.table$h12.hi, giant.table$h12.lo, giant.table$h12.hi)
giant.table$h24 <- ifelse(-giant.table$h24.lo>giant.table$h24.hi, giant.table$h24.lo, giant.table$h24.hi)
rownames(giant.table) <- giant.table$sig
selected.sigs <- c("SMID_BREAST_CANCER_BASAL_UP",
                   "SMID_BREAST_CANCER_BASAL_DN",
                   "REACTOME_CELL_CYCLE",
                   "ONDER_CDH1_TARGETS_2_DN",
                   "ONDER_CDH1_TARGETS_2_UP",
                   "ZWANG_CLASS_1_TRANSIENTLY_INDUCED_BY_EGF",
                   "DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP",
                   "DUTERTRE_ESTRADIOL_RESPONSE_24HR_DN",
                   "STK33_UP",
                   "KAUFFMANN_DNA_REPAIR_GENES",
                   "DANG_REGULATED_BY_MYC_DN",
                   "DANG_REGULATED_BY_MYC_UP"
                   )

heatmap.2(as.matrix(giant.table[giant.table$sig %in% selected.sigs, 10:13]),  
          breaks=seq(-20,20,2), trace="none", scale="none", col=colorRampPalette(c("blue", "white", "red")), margins=c(5,10), Colv=F, key=F)

heatmap.2(as.matrix(giant.table[, 10:13]),  
          breaks=seq(-20, 20, 2), trace="none", scale="none", col=colorRampPalette(c("blue", "white", "red")), margins=c(5,20), Colv=F, labRow="")





giant.hvw.sn <- s2n.rank(giant, giant.cls, list(c("h0n", "h6n", "h12n", "h24n"), c("w0n", "w6n", "w12n", "w24n")))
giant.hvw.table <- FES.table(giant.hvw.sn, sigs)
head(giant.hvw.table[order(giant.hvw.table$p.lo), ], 20)
head(giant.hvw.table[order(giant.hvw.table$p.hi), ], 20)
giant.hvw.table[grep("SMID", giant.hvw.table$sig), ]
giant.hvw.table[grep("CELL_CYCLE", giant.hvw.table$sig), ]
giant.hvw.table[grep("ONDER", giant.hvw.table$sig), ]
giant.hvw.table[grep("SREB", giant.hvw.table$sig), ]
giant.hvw.table[grep("DNA_REPAIR", giant.hvw.table$sig), ]
FES.plot(giant.hvw.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
FES.leading.edge(giant.hvw.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
FES.plot(giant.hvw.sn, sigs$gene[sigs$sig=="KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_UP"])
ES.plot(giant.hvw.sn, sigs$gene[sigs$sig=="KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_UP"])
write.csv(giant.hvw.table, "giant.hvw.csv")



giant.w6a.sn <- s2n.rank(giant, giant.cls, c("w6a", "w6n"))
giant.w6a.table <- FES.table(giant.w6a.sn, sigs)
head(giant.w6a.table[order(giant.w6a.table$p.lo), ], 20)
head(giant.w6a.table[order(giant.w6a.table$p.hi), ], 20)
giant.w6a.table[grep("SMID", giant.w6a.table$sig), ]

giant.h6r.sn <- s2n.rank(giant, giant.cls, c("h6r", "h6n"))
giant.h6r.table <- FES.table(giant.h6r.sn, sigs)
head(giant.h6r.table[order(giant.h6r.table$p.lo), ], 20)
head(giant.h6r.table[order(giant.h6r.table$p.hi), ], 20)

giant.w6r.sn <- s2n.rank(giant, giant.cls, c("w6r", "w6n"))
giant.w6r.table <- FES.table(giant.w6r.sn, sigs)
head(giant.w6r.table[order(giant.w6r.table$p.lo), ], 20)
head(giant.w6r.table[order(giant.w6r.table$p.hi), ], 20)

giant.w12r.sn <- s2n.rank(giant, giant.cls, c("w12r", "w12n"))
giant.w12r.table <- FES.table(giant.w12r.sn, sigs)
head(giant.w12r.table[order(giant.w12r.table$p.lo), ], 20)
head(giant.w12r.table[order(giant.w12r.table$p.hi), ], 20)

giant.w24r.sn <- s2n.rank(giant, giant.cls, c("w24r", "w24n"))
giant.w24r.table <- FES.table(giant.w24r.sn, sigs)
head(giant.w24r.table[order(giant.w24r.table$p.lo), ], 20)
head(giant.w24r.table[order(giant.w24r.table$p.hi), ], 20)

giant.w24a.sn <- s2n.rank(giant, giant.cls, c("w24a", "w24n"))
giant.w24a.table <- FES.table(giant.w24a.sn, sigs)
head(giant.w24a.table[order(giant.w24a.table$p.lo), ], 20)
head(giant.w24a.table[order(giant.w24a.table$p.hi), ], 20)

giant.6r.sn <- s2n.rank(giant, giant.cls, c("h6r", "w6r"))
giant.6r.table <- FES.table(giant.6r.sn, sigs)
head(giant.6r.table[order(giant.6r.table$p.lo), ], 20)
head(giant.6r.table[order(giant.6r.table$p.hi), ], 20)


# Normal GSEA
source("write.rnk.R")
write.sn.rnk(giant.0.sn, "giant.0.rnk")
write.sn.rnk(giant.6.sn, "giant.6.rnk")
write.sn.rnk(giant.12.sn, "giant.12.rnk")
write.sn.rnk(giant.24.sn, "giant.24.rnk")
all.0 <- GSEA.preranked.java.exec(input.rnk="giant.0.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "all.0")
all.6 <- GSEA.preranked.java.exec(input.rnk="giant.6.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "all.6")
all.12 <- GSEA.preranked.java.exec(input.rnk="giant.12.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "all.12")
all.24 <- GSEA.preranked.java.exec(input.rnk="giant.24.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "all.24")
gsea.all.data <- cbind(read.gsea.results(all.0), cat="c2.cgp", data="giant.0")
gsea.all.data <- rbind(gsea.all.data, cbind(read.gsea.results(all.6), cat="c2.cgp", data="giant.6"))
gsea.all.data <- rbind(gsea.all.data, cbind(read.gsea.results(all.12), cat="c2.cgp", data="giant.12"))
gsea.all.data <- rbind(gsea.all.data, cbind(read.gsea.results(all.24), cat="c2.cgp", data="giant.24"))

gsea.all.heatmap <- acast(gsea.all.data, sig ~ data, value.var="NES")
heatmap.2(gsea.all.heatmap, breaks=seq(-3,3,0.1), trace="none", scale="none", col=colorRampPalette(c("blue", "white", "red")), margins=c(5,10), Colv=F, key=F)

selected.sigs <- c("SMID_BREAST_CANCER_BASAL_UP",
                   "SMID_BREAST_CANCER_BASAL_DN",
                   "ONDER_CDH1_TARGETS_2_DN",
                   "ONDER_CDH1_TARGETS_2_UP",
                   "ZWANG_CLASS_1_TRANSIENTLY_INDUCED_BY_EGF",
                   "DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP",
                   "DUTERTRE_ESTRADIOL_RESPONSE_24HR_DN",
                   "KAUFFMANN_DNA_REPAIR_GENES",
                   "DANG_REGULATED_BY_MYC_DN",
                   "DANG_REGULATED_BY_MYC_UP"
)

png("gsea.heatmap.test.png", 1500, 1000)
heatmap.2(gsea.all.heatmap[selected.sigs, ], breaks=seq(-3,3,0.1), trace="none", scale="none", 
          col=colorRampPalette(c("blue", "white", "red")), margins=c(10,40), Colv=F, dendrogram="none", key=F,
          cexRow=1.5, cexCol=1.5)
dev.off()

ES.plot(giant.0.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
ES.plot(giant.6.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
ES.plot(giant.12.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
ES.plot(giant.24.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])

png("cell.cycle.heatmap.png", 1500, 1000)
heatmap.2(as.matrix(giant.scaled[giant.scaled$NAME %in% sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"], c(3:5, 27:29, 9:11, 18:20, 33:35, 57:60, 39:41, 48:50)]), 
          breaks=seq(-3,3,0.1), trace="none", col=colorRampPalette(c("blue", "white", "red")), Colv=F, labRow="", dendrogram="none", key=F,
          cexCol=1.5)
dev.off()

png("global.heatmap.png", 1500, 1000)
heatmap.2(as.matrix(giant.scaled[, c(3:5, 27:29, 9:11, 18:20, 33:35, 57:60, 39:41, 48:50)]), 
          breaks=seq(-3,3,0.1), trace="none", col=colorRampPalette(c("blue", "white", "red")),
          cexCol=1.5)
dev.off()

write.sn.rnk(tcga.bvn.sn, "tcga.bvn.rnk")
c1.bvn  <- GSEA.preranked.java.exec(input.rnk="tcga.bvn.rnk", gs.db="c1.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c1.bvn")
c2.cp.bvn <- GSEA.preranked.java.exec(input.rnk="tcga.bvn.rnk", gs.db="c2.cp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cp.bvn")
c2.cgp.bvn <- GSEA.preranked.java.exec(input.rnk="tcga.bvn.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cgp.bvn")
c3.bvn <- GSEA.preranked.java.exec(input.rnk="tcga.bvn.rnk", gs.db="c3.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c3.bvn")
c5.bvn <- GSEA.preranked.java.exec(input.rnk="tcga.bvn.rnk", gs.db="c5.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c5.bvn")
c6.bvn <- GSEA.preranked.java.exec(input.rnk="tcga.bvn.rnk", gs.db="c6.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c6.bvn")

write.sn.rnk(tcga.bvl.sn, "tcga.bvl.rnk")
c1.bvl  <- GSEA.preranked.java.exec(input.rnk="tcga.bvl.rnk", gs.db="c1.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c1.bvl")
c2.cp.bvl <- GSEA.preranked.java.exec(input.rnk="tcga.bvl.rnk", gs.db="c2.cp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cp.bvl")
c2.cgp.bvl <- GSEA.preranked.java.exec(input.rnk="tcga.bvl.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cgp.bvl")
c3.bvl <- GSEA.preranked.java.exec(input.rnk="tcga.bvl.rnk", gs.db="c3.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c3.bvl")
c5.bvl <- GSEA.preranked.java.exec(input.rnk="tcga.bvl.rnk", gs.db="c5.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c5.bvl")
c6.bvl <- GSEA.preranked.java.exec(input.rnk="tcga.bvl.rnk", gs.db="c6.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c6.bvl")

write.sn.rnk(mtor0.rna.sn.scaled, "mtor0.rna.rnk")
c1.rna  <- GSEA.preranked.java.exec(input.rnk="mtor0.rna.rnk", gs.db="c1.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c1.rna")
c2.cp.rna <- GSEA.preranked.java.exec(input.rnk="mtor0.rna.rnk", gs.db="c2.cp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cp.rna")
c2.cgp.rna <- GSEA.preranked.java.exec(input.rnk="mtor0.rna.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cgp.rna")
c3.rna <- GSEA.preranked.java.exec(input.rnk="mtor0.rna.rnk", gs.db="c3.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c3.rna")
c5.rna <- GSEA.preranked.java.exec(input.rnk="mtor0.rna.rnk", gs.db="c5.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c5.rna")
c6.rna <- GSEA.preranked.java.exec(input.rnk="mtor0.rna.rnk", gs.db="c6.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c6.rna")

write.sn.rnk(mtor0.silac.sn, "mtor0.prot.rnk")
c1.prot  <- GSEA.preranked.java.exec(input.rnk="mtor0.prot.rnk", gs.db="c1.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c1.prot")
c2.cp.prot <- GSEA.preranked.java.exec(input.rnk="mtor0.prot.rnk", gs.db="c2.cp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cp.prot")
c2.cgp.prot <- GSEA.preranked.java.exec(input.rnk="mtor0.prot.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cgp.prot")
c3.prot <- GSEA.preranked.java.exec(input.rnk="mtor0.prot.rnk", gs.db="c3.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c3.prot")
c5.prot <- GSEA.preranked.java.exec(input.rnk="mtor0.prot.rnk", gs.db="c5.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c5.prot")
c6.prot <- GSEA.preranked.java.exec(input.rnk="mtor0.prot.rnk", gs.db="c6.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c6.prot")

write.sn.rnk(tcga.pik3ca.sn, "tcga.b1047.rnk")
c1.b1047  <- GSEA.preranked.java.exec(input.rnk="tcga.b1047.rnk", gs.db="c1.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c1.b1047")
c2.cp.b1047 <- GSEA.preranked.java.exec(input.rnk="tcga.b1047.rnk", gs.db="c2.cp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cp.b1047")
c2.cgp.b1047 <- GSEA.preranked.java.exec(input.rnk="tcga.b1047.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cgp.b1047")
c3.b1047 <- GSEA.preranked.java.exec(input.rnk="tcga.b1047.rnk", gs.db="c3.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c3.b1047")
c5.b1047 <- GSEA.preranked.java.exec(input.rnk="tcga.b1047.rnk", gs.db="c5.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c5.b1047")
c6.b1047 <- GSEA.preranked.java.exec(input.rnk="tcga.b1047.rnk", gs.db="c6.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c6.b1047")

tcga.pik3ca <- tcga.add.maf(tcga.txt, tcga.maf, "PIK3CA")
tcga.pik3ca.only.cls <- tcga.cls(tcga.breast, tcga.pik3ca, tcga.pik3ca$Sample[tcga.pik3ca$PIK3CA], tcga.pik3ca$Sample[tcga.pik3ca$PAM50 == "Normal"])
tcga.pik3ca.only.sn <- s2n.rank(tcga.breast, tcga.pik3ca.only.cls)
write.sn.rnk(tcga.pik3ca.only.sn, "tcga.1047.rnk")
c1.1047  <- GSEA.preranked.java.exec(input.rnk="tcga.1047.rnk", gs.db="c1.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c1.1047")
c2.cp.1047 <- GSEA.preranked.java.exec(input.rnk="tcga.1047.rnk", gs.db="c2.cp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cp.1047")
c2.cgp.1047 <- GSEA.preranked.java.exec(input.rnk="tcga.1047.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cgp.1047")
c3.1047 <- GSEA.preranked.java.exec(input.rnk="tcga.1047.rnk", gs.db="c3.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c3.1047")
c5.1047 <- GSEA.preranked.java.exec(input.rnk="tcga.1047.rnk", gs.db="c5.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c5.1047")
c6.1047 <- GSEA.preranked.java.exec(input.rnk="tcga.1047.rnk", gs.db="c6.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c6.1047")


write.sn.rnk(giant.hvw.sn, "giant.hvw.rnk")
c1.hvw  <- GSEA.preranked.java.exec(input.rnk="giant.hvw.rnk", gs.db="c1.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c1.hvw")
c2.cp.hvw  <- GSEA.preranked.java.exec(input.rnk="giant.hvw.rnk", gs.db="c2.cp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cp.hvw")
c2.cgp.hvw  <- GSEA.preranked.java.exec(input.rnk="giant.hvw.rnk", gs.db="c2.cgp.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c2.cgp.hvw")
c3.hvw  <- GSEA.preranked.java.exec(input.rnk="giant.hvw.rnk", gs.db="c3.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c3.hvw")
c5.hvw  <- GSEA.preranked.java.exec(input.rnk="giant.hvw.rnk", gs.db="c5.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c5.hvw")
c6.hvw  <- GSEA.preranked.java.exec(input.rnk="giant.hvw.rnk", gs.db="c6.all.v4.0.symbols.gmt",  output.directory = "gsea",  gs.size.threshold.max=1500, doc.string = "c6.hvw")



gsea.data <- cbind(read.gsea.results(c1.1047), cat="c1", data="1047")
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cp.1047), cat="c2.cp", data="1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cgp.1047), cat="c2.cgp", data="1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c3.1047), cat="c3", data="1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.1047), cat="c5", data="1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.1047), cat="c6", data="1047"))

gsea.data <- rbind(gsea.data,cbind(read.gsea.results(c1.b1047), cat="c1", data="b1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cp.b1047), cat="c2.cp", data="b1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cgp.b1047), cat="c2.cgp", data="b1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c3.b1047), cat="c3", data="b1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.b1047), cat="c5", data="b1047"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.b1047), cat="c6", data="b1047"))

gsea.data <- rbind(gsea.data,cbind(read.gsea.results(c1.prot), cat="c1", data="prot"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cp.prot), cat="c2.cp", data="prot"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cgp.prot), cat="c2.cgp", data="prot"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c3.prot), cat="c3", data="prot"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.prot), cat="c5", data="prot"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.prot), cat="c6", data="prot"))

gsea.data <- rbind(gsea.data,cbind(read.gsea.results(c1.rna), cat="c1", data="rna"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cp.rna), cat="c2.cp", data="rna"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cgp.rna), cat="c2.cgp", data="rna"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c3.rna), cat="c3", data="rna"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.rna), cat="c5", data="rna"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.rna), cat="c6", data="rna"))

gsea.data <- rbind(gsea.data,cbind(read.gsea.results(c1.bvl), cat="c1", data="bvl"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cp.bvl), cat="c2.cp", data="bvl"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cgp.bvl), cat="c2.cgp", data="bvl"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c3.bvl), cat="c3", data="bvl"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.bvl), cat="c5", data="bvl"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.bvl), cat="c6", data="bvl"))

gsea.data <- rbind(gsea.data,cbind(read.gsea.results(c1.bvn), cat="c1", data="bvn"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cp.bvn), cat="c2.cp", data="bvn"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c2.cgp.bvn), cat="c2.cgp", data="bvn"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c3.bvn), cat="c3", data="bvn"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.bvn), cat="c5", data="bvn"))
gsea.data <- rbind(gsea.data, cbind(read.gsea.results(c5.bvn), cat="c6", data="bvn"))

# tile view
selected.sigs <- c(#"DANG_REGULATED_BY_MYC_DN",
                   "DANG_REGULATED_BY_MYC_UP",
                   #"DANG_MYC_TARGETS_UP",
                   "CELL_CYCLE_PHASE",
                   "KAUFFMANN_DNA_REPAIR_GENES",
                   #"ONDER_CDH1_TARGETS_3_DN",
                   #"ONDER_CDH1_TARGETS_2_DN",
                   #"ONDER_CDH1_TARGETS_1_DN",
                   #"ONDER_CDH1_TARGETS_3_UP",
                   "ONDER_CDH1_TARGETS_2_UP",
                   #"ONDER_CDH1_TARGETS_1_UP",
                   "SMID_BREAST_CANCER_BASAL_UP"
                   #"SMID_BREAST_CANCER_BASAL_DN"
                   
)
selected.data <- c("rna", "bvn", "b1047")
temp.data <- gsea.data[gsea.data$sig %in% selected.sigs & gsea.data$data %in% selected.data, ]
data2 <- data.frame(data=c("bvn", "b1047", "rna"), data2=c("Basal vs Normal", "Basal w PIK3CA mutation vs Normal", "MCF-10A H1047R vs WT"))
temp.data <- merge(temp.data, data2, by="data")
sig2 <- data.frame(sig=c("SMID_BREAST_CANCER_BASAL_UP", "DANG_REGULATED_BY_MYC_UP", "CELL_CYCLE_PHASE", "KAUFFMANN_DNA_REPAIR_GENES", 
                     "ONDER_CDH1_TARGETS_2_UP"), 
                   sig2=c("Basal Breast Cancer", "Myc Upregulated", "Cell cycle", "DNA Repair", "Beta Catenin"))
temp.data <- merge(temp.data, sig2, by="sig")

png("table.png", 2000, 1500, res=300)
ggplot(temp.data, aes(y=sig2, x=data2, fill=-log10(p)*sign(NES))) +
  geom_tile() + scale_fill_gradient2(high="#BB2020", low="#2020BB", name="Up/Down") +
  theme_bw()  + theme(text=element_text(size=8)) + ylab("") + xlab("")
dev.off()

png("table-NES.png", 2000, 1500, res=300)
ggplot(temp.data, aes(y=sig2, x=data2, fill=NES)) +
  geom_tile() + scale_fill_gradient2(high="#BB2020", low="#2020BB", name="NES") +
  theme_bw()  + theme(text=element_text(size=8)) + ylab("") + xlab("")
dev.off()

# figures
png("charafe.tcga.bvl.png", 4,5,units="in",res=300)
ES.plot(tcga.bvl.sn, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN"])
dev.off()

png("charafe.tcga.bvn.png", 4,5,units="in",res=300)
ES.plot(tcga.bvn.sn, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN"])
dev.off()

png("smid.tcga.bvl.png", 4,5,units="in",res=300)
FES.plot(tcga.bvl.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()

png("smid.tcga.bvn.png", 4,5,units="in",res=300)
FES.plot(tcga.bvn.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()

png("smid.rna.png", 4,5,units="in",res=300)
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()
FES.leading.edge(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])


png("smid.silac.png", 4,5,units="in",res=300)
FES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()

png("smid.b1047.png", 4,5,units="in",res=300)
FES.plot(tcga.pik3ca.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()

png("onder.silac.png", 4,5,units="in",res=300)
ES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
dev.off()

png("onder.rna.png", 4,5,units="in",res=300)
ES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
dev.off()

png("onder.bvl.png", 4,5,units="in",res=300)
ES.plot(tcga.bvl.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
dev.off()

png("onder.bvn.png", 4,5,units="in",res=300)
ES.plot(tcga.bvn.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
dev.off()

png("onder.b1047.png", 4,5,units="in",res=300)
ES.plot(tcga.pik3ca.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
dev.off()

png("cell.cycle.silac.png", 4,5,units="in",res=300)
ES.plot(mtor0.silac.sn, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])
dev.off()

png("cell.cycle.rna.png", 4,5,units="in",res=300)
ES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])
dev.off()

png("cell.cycle.bvl.png", 4,5,units="in",res=300)
ES.plot(tcga.bvl.sn, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])
dev.off()

png("cell.cycle.bvn.png", 4,5,units="in",res=300)
ES.plot(tcga.bvn.sn, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])
dev.off()

png("cell.cycle.b1047.png", 4,5,units="in",res=300)
ES.plot(tcga.pik3ca.sn, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])
dev.off()

png("fes.smid.rna.png", 4,5,units="in",res=300)
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()

png("fes.smid.bvl.png", 4,5,units="in",res=300)
FES.plot(tcga.bvl.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()

png("fes.smid.bvn.png", 4,5,units="in",res=300)
FES.plot(tcga.bvn.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
dev.off()

write.csv(mtor0.rna.table.scaled, "mtor0.rna.csv")
write.csv(mtor0.silac.table, "mtor0.silac.csv")
write.csv(tcga.bvl.table, "tcga.bvl.csv")
write.csv(tcga.bvn.table, "tcga.bvn.csv")
write.csv(tcga.pik3ca.table, "tcga.b1047.csv")

FES.plot(tcga.bvl.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="ONDER_CDH1_TARGETS_2_DN"])
FES.plot(tcga.bvn.sn, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])

# PIK3CA mutant breast vs not
tcga.pik3ca <- tcga.add.maf(tcga.txt, tcga.maf, "PIK3CA")
tcga.pik3ca.cls <- tcga.cls(tcga.breast, tcga.pik3ca, tcga.pik3ca$Sample[tcga.pik3ca$PIK3CA & tcga.pik3ca$PAM50 == "Basal"], tcga.pik3ca$Sample[tcga.pik3ca$PAM50 == "Normal"])
tcga.pik3ca.sn <- s2n.rank(tcga.breast, tcga.pik3ca.cls)
tcga.pik3ca.table <- FES.table(tcga.pik3ca.sn, sigs)
head(tcga.pik3ca.table[order(tcga.pik3ca.table$p.hi),], 20)
head(tcga.pik3ca.table[order(tcga.pik3ca.table$p.lo),], 20)
FES.plot(tcga.pik3ca.sn, sigs$gene[sigs$sig=="CHANG_CYCLING_GENES"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="CHANG_CYCLING_GENES"])
FES.plot(tcga.pik3ca.sn, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])
FES.plot(mtor0.rna.sn.scaled, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])

tcga.clust <- hclust(dist(as.matrix(tcga.breast[, c(-1,-2)])))
# random sets test
genes <- unique(sigs$gene)
set.sizes <- ddply(sigs, .(sig), nrow)
set.seed(123)
random.sets <- lapply(1:10000, function(x) sample.int(length(genes), set.sizes$V1[round(runif(1, 1, nrow(set.sizes)))]))
random.sigs <- foreach(i = 1:length(random.sets), .combine=rbind) %do% {
  data.frame(sig=paste0("rand",i), gene=genes[random.sets[[i]]])
}
mtor0.rna.rand.table <- FES.table(mtor0.rna.sn, random.sigs)
hist(-log10(mtor0.rna.rand.table$p.hi))
hist(-log10(mtor0.rna.rand.table$p.lo))
hist(mtor0.rna.rand.table$q.hi)
hist(mtor0.rna.rand.table$q.lo)
head(mtor0.rna.rand.table[order(mtor0.rna.rand.table$p.hi),], 20)
head(mtor0.rna.rand.table[order(mtor0.rna.rand.table$p.lo),], 20)
FES.plot(mtor0.rna.sn, random.sigs$gene[random.sigs$sig=="rand8918"])
