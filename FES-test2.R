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
ggplot(p4936.sn.dox, aes(x=rank, y=sn)) + geom_point()
mtor0.silac.table <- FES.table(mtor0.silac.sn, sigs)
head(mtor0.silac.table[order(mtor0.silac.table$p.hi),], 20)
head(mtor0.silac.table[order(mtor0.silac.table$p.lo),], 20)
sum(mtor0.silac.table$q.hi<0.01)
sum(mtor0.silac.table$q.lo<0.01)
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
mtor0.rna.cls <- read.cls("datasets//hvw.cls")
mtor0.rna.sn <- s2n.rank(mtor0.rna, mtor0.rna.cls)
mtor0.rna.table <- FES.table(mtor0.rna.sn, sigs)
head(mtor0.rna.table[order(mtor0.rna.table$p.hi),], 20)
head(mtor0.rna.table[order(mtor0.rna.table$p.lo),], 20)
sum(mtor0.rna.table$q.hi<0.01)
sum(mtor0.rna.table$q.lo<0.01)
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="REACTOME_CELL_CYCLE"])
foo <- FES.leading.edge(mtor0.rna.sn, sigs$gene[sigs$sig=="HUPER_BREAST_BASAL_VS_LUMINAL_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="HUPER_BREAST_BASAL_VS_LUMINAL_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="HUPER_BREAST_BASAL_VS_LUMINAL_DN"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="SMID_BREAST_CANCER_BASAL_DN"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="MYC_UP.V1_DN"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="MYC_UP.V1_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="ODONNELL_TARGETS_OF_MYC_AND_TFRC_UP"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="ODONNELL_TARGETS_OF_MYC_AND_TFRC_DN"])
FES.plot(p4936.sn.dox, sigs$gene[sigs$sig=="KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_DN"])

# TCGA breast data
tcga.breast <- read.tcga("BRCA.exp.547.med.txt")
tcga.cls <- read.tcga.cls("BRCA.547.PAM50.SigClust.Subtypes.txt", "PAM50", tcga.breast)
tcga.maf <- read.maf("genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf")
tcga.txt <- read.tcga.txt("BRCA.547.PAM50.SigClust.Subtypes.txt")
with(tcga.maf, tcga.maf[Hugo_Symbol=="PIK3CA" & amino_acid_change_WU=="p.H1047R", ])
# Basal will be positive(high) Luminal will be negative(low)
tcga.sn <- s2n.rank(tcga.breast, tcga.cls, c("Basal", "LumA"))
tcga.bvl.table <- FES.table(tcga.sn, sigs)
head(tcga.bvl.table[order(tcga.table$p.hi),], 20)
head(tcga.bvl.table[order(tcga.table$p.lo),], 20)

# Basal will be positive(high) normal will be negative(low)
tcga.sn <- s2n.rank(tcga.breast, tcga.cls, c("Basal", "Normal"))
tcga.bvn.table <- FES.table(tcga.sn, sigs)
head(tcga.bvl.table[order(tcga.table$p.hi),], 20)
head(tcga.bvl.table[order(tcga.table$p.lo),], 20)
FES.plot(tcga.sn, sigs$gene[sigs$sig=="PUJANA_BRCA1_PCC_NETWORK"])
FES.plot(mtor0.rna.sn, sigs$gene[sigs$sig=="PUJANA_BRCA1_PCC_NETWORK"])

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
