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
require(plyr)
require(doMC)
require(ggplot2)
require(grid)
require(gridExtra)
require(reshape2)
# Determine signal to noise scoring for each line in the data table
s2n.score <- function(dataA, dataB) {
  mu1 <- apply(dataA, 1, mean)
  mu2 <- apply(dataB, 1, mean)
  sigma1 <- apply(dataA, 1, sd)
  sigma2 <- apply(dataB, 1, sd)
  (mu1-mu2)/(sigma1+sigma2)
}
filter.cpm <- function(gct, cls, min.cpm=1, contrast=NULL) {
  if (is.null(contrast)){
    contrast=list(levels(cls)[1], levels(cls)[2])
  }
  len1 <- sum(cls %in% contrast[[1]])
  len2 <- sum(cls %in% contrast[[2]])
  selection1 <- rowSums(gct[, which(cls %in% contrast[[1]])+2]>min.cpm)
  selection2 <- rowSums(gct[, which(cls %in% contrast[[2]])+2]>min.cpm)
  gct[(selection1==len1) & (selection2==len2),]
}

s2n.rank <- function(gct, cls, contrast=NULL) {
  if (is.null(contrast)){
    contrast=list(levels(cls)[1], levels(cls)[2])
  }
  sn.vector <- s2n.score(gct[, which(cls %in% contrast[[1]])+2], gct[, which(cls %in% contrast[[2]])+2])
  sn.vector <- data.frame(sn=sn.vector)
  rownames(sn.vector) <- gct[, 1]
  sn.vector <- na.omit(sn.vector)
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector <- sn.vector[order(-sn.vector$sn), ]
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector
}

z.scale <- function(gct, cls) {
  exp.matrix <- as.matrix(gct[, c(-1, -2)])
  exp.matrix <- log2(exp.matrix)
  # row medians
  median.vector <- apply(exp.matrix, 1, median)
  # sd by sample group and z-scale
  scaled.matrix <- foreach(sample.group=iter(unique(cls)), .combine=cbind) %do% {
    sd.vector <- apply(exp.matrix[, which(cls == sample.group)], 1, sd)
    (exp.matrix[, which(cls == sample.group)] - median.vector)/sd.vector
#     sd.vector <- data.frame(sd.vector)
#     colnames(sd.vector) <- sample.group
#     sd.vector
  }
  data.frame(gct[, c(1,2)], scaled.matrix)
  
}

z.score <- function(gct, cls) {
  exp.matrix <- as.matrix(gct[, c(-1, -2)])
  exp.matrix <- log2(exp.matrix)
  mean.int <- apply(exp.matrix, 2, mean)
  sd.int <- apply(exp.matrix, 2, sd)
  z.score.mat <- foreach(i=1:ncol(exp.matrix), .combine=cbind) %do% {
    exp.matrix[, i] <- (exp.matrix[, i] - mean.int[i])/sd.int[i]
  }
  colnames(z.score.mat) <- colnames(exp.matrix)
  data.frame(gct[, c(1,2)], z.score.mat)
  
}
z.ratio <- function(dataA, dataB) {
  mu1 <- apply(dataA, 1, mean)
  mu2 <- apply(dataB, 1, mean)
  sigma1 <- apply(dataA, 1, sd)
  sigma2 <- apply(dataB, 1, sd)
  n1 <- ncol(dataA)
  n2 <- ncol(dataB)
  (mu1-mu2)/sqrt((sigma1^2/n1)+(sigma2^2/n2))
}
z.rank <- function(gct, cls, contrast=NULL) {
  if (is.null(contrast)){
    contrast=list(levels(cls)[1], levels(cls)[2])
  }
  sn.vector <- z.ratio(gct[, which(cls %in% contrast[[1]])+2], gct[, which(cls %in% contrast[[2]])+2])
  sn.vector <- data.frame(sn=sn.vector)
  rownames(sn.vector) <- gct[, 1]
  sn.vector <- na.omit(sn.vector)
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector <- sn.vector[order(-sn.vector$sn), ]
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector
}
  
noise.rank <- function(gct, cls, contrast=NULL) {
  if (is.null(contrast)){
    contrast=list(levels(cls)[1])
  }
  sd.vector <- apply(gct[, which(cls %in% contrast[[1]])+2], 1, sd)
  sn.vector <- data.frame(sn=sd.vector)
  rownames(sn.vector) <- gct[, 1]
  sn.vector <- na.omit(sn.vector)
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector <- sn.vector[order(-sn.vector$sn), ]
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector
}  


s2n.direct.rank <- function(gct, high.samples, low.samples) {
  sn.vector <- s2n.score(gct[, which(cls %in% contrast[[1]])+2], gct[, which(cls %in% contrast[[2]])+2])
  sn.vector <- data.frame(sn=sn.vector)
  rownames(sn.vector) <- gct[, 1]
  sn.vector <- na.omit(sn.vector)
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector <- sn.vector[order(-sn.vector$sn), ]
  sn.vector$rank <- 1:nrow(sn.vector)
  sn.vector
}

FES <- function(sn.table, geneset, hits.only=TRUE) {
  genes <- rownames(sn.table)
  overlap <- intersect(genes, geneset)
  hit <- genes %in% overlap
  # preallocate
  lower <- vector("numeric", nrow(sn.table))
  higher <- vector("numeric", nrow(sn.table))
  hit.count <- 0
  miss.count <- 0
  for(index in 1:nrow(sn.table)) {
    if(hit[index]) {
      hit.count <- hit.count + 1
      lower[index] <- phyper(hit.count, length(overlap), length(genes)-length(overlap), index) #* (length(overlap) - hit.count + 1)
      higher[index] <- phyper(length(overlap)-hit.count, length(overlap), length(genes)-length(overlap), length(genes)-index) #*  (hit.count + 1)
    } else {
      miss.count <- miss.count + 1
      if (hits.only) {
        lower[index] <- 1
        higher[index] <- 1
      } else {
        lower[index] <- phyper(hit.count, length(overlap), length(genes)-length(overlap), index) #* (length(overlap) - hit.count + 1)
        higher[index] <- phyper(length(overlap)-hit.count, length(overlap), length(genes)-length(overlap), length(genes)-index) #*  (hit.count + 1)
      }
    }
  }
  lower[lower>1] <- 1
  higher[higher>1] <- 1
  data.frame(sn.table, p.lower=lower, p.higher=higher, hit=hit)
}
FES.summary <- function(sn.table, geneset) {
  fes.data <- FES(sn.table, geneset)
  # p.value correction for multiple testing
  fes.data$p.lower <- fes.data$p.lower
  with(fes.data, data.frame(overlap=sum(hit),
                            max.lo=which.min(p.lower), 
                            sn.lo=sn[which.min(p.lower)], 
                            p.lo=min(p.lower), 
                            max.hi=which.min(p.higher), 
                            sn.hi=sn[which.min(p.higher)], 
                            p.hi=min(p.higher)))
}
FES.table <- function(sn.table, genesets) {
  # this factor is to drop unused levels from the sig list which causes an issue with ddply
  genesets$sig <- factor(genesets$sig)
  out.table <- ddply(genesets, .(sig), .fun=function(df) FES.summary(sn.table, df$gene), .parallel=T)
  out.table$q.lo <- p.adjust(out.table$p.lo, method="BH")
  out.table$q.hi <- p.adjust(out.table$p.hi, method="BH")
  out.table[, c(1:5, 9, 6:8, 10)]
}
FES.plot <- function(sn.table, geneset) {
  
  plot.data <- FES(sn.table, geneset, hits.only=F) 
  plot.data$percentile <- as.numeric(cut(plot.data$sn, quantile(plot.data$sn, probs=seq(0, 1, 1/9)), include.lowest=T))
  plot.data$hit <- ifelse(plot.data$hit, 1, 0)

  
  plot.data2 <- melt(plot.data, id.vars="rank")
  
  # need to combine the two p value data tracks into one
  plot.data2$variable <- as.character(plot.data2$variable)
  plot.data2$hilo <- NA
  plot.data2$hilo[plot.data2$variable == "p.lower"] <- "lower"
  plot.data2$hilo[plot.data2$variable == "p.higher"] <- "higher"
  plot.data2$hilo <- factor(plot.data2$hilo)
  plot.data2$variable[plot.data2$variable %in% c("p.higher", "p.lower")] <- "p"
  
  peak.pos <- plot.data2$rank[which.min(plot.data2$value[plot.data2$variable == "p"])]
  
  # fix the facet order
  plot.data2$variable <- factor(plot.data2$variable, levels=c("p", "hit", "percentile", "sn"))
  
  p <- ggplot(plot.data2, aes(x=rank)) + 
    facet_wrap(~ variable, ncol=1, scales="free_y") +
    geom_line(subset=.(variable == "p"), aes(y=-log10(value), color=hilo)) + scale_color_manual(values=c("red","blue")) +
    
    geom_point(subset=.(variable == "sn"), aes(y=value)) +
    geom_segment(subset=.(variable == "hit"), aes(x=rank, xend=rank, y=value, yend=0)) +
    geom_tile(subset=.(variable == "percentile"), aes(y=1, fill=value, width=2)) + scale_fill_gradient2(midpoint=5, low="blue", high="red") +
    #geom_rect(aes(xmin=index.min, xmax=index.max, color=index), ymin=0, ymax=1, data=percentile.data) +
    #geom_segment(subset=.(variable == "percentile"), aes(y=1, yend=0, x=rank, xend=rank, color=value)) + scale_color_gradient2(midpoint=5, low="blue", high="red") +
    geom_vline(xintercept=peak.pos, color="darkgreen") +
    theme_bw() +
    theme(strip.text=element_blank(), 
          legend.position="none", 
          strip.background = element_blank(),
          plot.margin = unit( c(0,0,0,0) , units = "lines"),
          panel.margin= unit(c(0,0,0,0), units="null"),
          panel.border=element_blank()) +
    ylab("")
   
  # hacky things to modify the facets to different heights, remove the axis that aren't needed, etc.
  require(grid)
  require(gtable)
  gt <- ggplot_gtable(ggplot_build(p))
  
  gt$heights[[4]] <- unit(3, "null") 
  gt$heights[[16]] <- unit(2, "null")
  gt <- gtable_add_grob(gt, rectGrob(gp=gpar(fill="white", col="white")), 7, 3, 15, name="blank1")
  gt <- gtable_add_grob(gt, textGrob(expression(-log[10]*p), rot=90, gp=gpar(fontsize=11)), 4, 2, name="blank1")
  gt <- gtable_add_grob(gt, textGrob("Signal to Noise", rot=90, gp=gpar(fontsize=11)), 16, 2, name="blank1")
  grid.draw(gt) 
}

FES.plot2 <- function(sn.table, geneset) {
  
  plot.data <- FES(sn.table, geneset, hits.only=F) 
  
  plot.data$percentile <- cut(plot.data$sn, quantile(plot.data$sn, probs=seq(0, 1, 0.11111)))
  #ggplot(plot.data, aes(x=rank)) + geom_line(aes(y=-log10(p.lower)), color="blue") + geom_line(aes(y=-log10(p.higher)), color="red")
  
  grid.newpage()
  sn.plot <- ggplot(plot.data, aes(x=rank, y=sn)) + geom_point() + theme_bw() +
    theme(plot.margin=unit(c(0, 1, 0, 0), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  ES.plot <- ggplot(plot.data, aes(x=rank, )) + geom_line(aes(y=-log10(p.lower)), color="blue") +
    geom_line(aes(y=-log10(p.higher)), color="red") +
    theme_bw() +
    theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(plot.margin=unit(c(0, 1, 0, 0), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  ticks <- ggplot(plot.data, aes(x=rank, xend=rank, y=ifelse(hit, 0.1, 0), yend=0)) +geom_segment() + theme_bw() +
    theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.line=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(plot.margin=unit(c(0, 1, 0, 1.1), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  gradient <- ggplot(plot.data, aes(x=rank))  +theme_bw() +
    geom_tile(aes(y=0, fill=as.numeric(percentile))) + scale_fill_gradient2(midpoint=5, low="blue", high="red") + theme(legend.position="none") +
    theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.line=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(plot.margin=unit(c(0, 1, 0, 1.1), "lines"), panel.border=element_rect(size=0)) +
    theme(panel.grid=element_blank())
  
  # hacky things to modify the facets to different heights

  
  grid.arrange(ES.plot, ticks, gradient, sn.plot, ncol=1, nrow=4, heights=c(5, 1, 1, 3))
  
  
}

FES.leading.edge <- function(sn.table, geneset) {
  results <- FES(sn.table, geneset)
  peak.p.lo <- results$p.lower[which.min(results$p.lower)]
  peak.p.hi <- results$p.higher[which.min(results$p.higher)]
  hilo <- peak.p.hi < peak.p.lo # hilo will be true if p.hi is less than p.lo aka hi is more significant
  if (hilo) {
    peak.hi <- results$rank[which.min(results$p.higher)]
    leading.edge <- results[results$hit & (results$rank <= peak.hi), ]
  } else {
    peak.lo <- results$rank[which.min(results$p.lower)]
    leading.edge <- results[results$hit & (results$rank >= peak.lo), ]
  }
  leading.edge
}

sig.overlap <- function(sig1, sig2, sig.table) {
  genes1 <- sig.table$gene[sig.table$sig==sig1]
  genes2 <- sig.table$gene[sig.table$sig==sig2]
  length(intersect(genes1, genes2)) / length(union(genes1, genes2))
}
sig.overlap.plot <- function(in.sigs, sig.table) {
  out.matrix <- matrix(data=vector(mode="numeric", length=length(in.sigs)^2), nrow=length(in.sigs), ncol=length(in.sigs))
  for(i in 1:length(in.sigs)) {
    for(j in 1:length(in.sigs)) {
      out.matrix[i, j] <- sig.overlap(in.sigs[i], in.sigs[j], sig.table)
    }
  }
  heatmap(out.matrix)
}

FES.heatmap <- function(gct, sigs) {
  heatmap.data <- foreach(i=3:ncol(gct), .combine=cbind) %do% {
    # rank the selected column
    exp.sample <- data.frame(sn=gct[, i], row.names=gct[, 1])
    exp.sample$rank <- 1
    exp.sample <- exp.sample[order(-exp.sample$sn), ]
    exp.sample$rank <- 1:nrow(exp.sample)
    
    out.table <- FES.table(exp.sample, sigs)
    
    # decide tests and magnitudes
    out.vector <- vector(mode="numeric", length=nrow(out.table))
    out.vector <- apply(out.table[, c(5, 9)], 1, min)
    out.vector <- -log10(out.vector)
    out.vector <- out.vector * ifelse(out.table$p.lo < out.table$p.hi, -1, 1)
    out.vector <- data.frame(out.vector)
    colnames(out.vector) <- colnames(gct)[i]
    out.vector
  }
  heatmap.data <- as.matrix(heatmap.data)
  rownames(heatmap.data) <- levels(factor(sigs$sig))
  heatmap.data
  #heatmap(heatmap.data)
}