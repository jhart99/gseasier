GSEA.java.exec <- function(
  input.ds,           # Input gene expression Affy dataset file in RES or GCT format
  input.cls,           # Input class vector (phenotype) file in CLS format
  gs.db,         # Gene set database in GMT format
  output.directory = "",        # Directory where to store output and results (default: "")
  doc.string            = "GSEA.analysis",   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
  non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
  reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
  nperm                 = 1000,            # Number of random permutations (default: 1000)
  weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
  nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
  fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
  fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
  topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
  adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
  gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
  gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
  reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
  preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
  random.seed           = 3338,            # Random number generator seed. (default: 123456)
  perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
  fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
  replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
  save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
  OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
  use.fast.enrichment.routine = T,          # Use faster routine to compute enrichment for random permutations (default: T)
  metric = "Signal2Noise",
  contrast = "",
  silent = T
){
  input.cls <- if (contrast !="") paste0(input.cls, "#", contrast)
  system(paste("java -cp ~/gsea/gsea2-2.0.13.jar ",
    "-Xmx4000m xtools.gsea.Gsea ",
    "-res ", input.ds,  " ",
    "-cls ", input.cls, " ",
    "-gmx ", gs.db, " ", 
    "-collapse false ",
    "-norm meandiv ",
    "-nperm 1000 ",
    "-permute gene_set ",
    "-rnd_type no_balance ",
    "-scoring_scheme weighted ",
    "-rpt_label ", doc.string, " ",
    "-metric ", metric, " ",
    "-num 100 ", 
    "-plot_top_x ", topgs, " ",
    "-rnd_seed ", random.seed, " ",
    "-save_rnd_lists false ",
    "-set_max ", gs.size.threshold.max, " ",
    "-set_min ", gs.size.threshold.min, " ", 
    "-zip_report false ",
    "-out ", output.directory, " ",
    "-gui false"), intern=silent)
  report.dir <- list.files(output.directory)
  file.index <- max(grep(paste0("^",doc.string), report.dir))
  # Return the output files
  list.files(paste0(output.directory, "/", report.dir[file.index]), "gsea_report.*html", full.names=T)
}
GSEA.preranked.java.exec <- function(
  input.rnk,           # Input gene expression Affy dataset file in RES or GCT format
  gs.db,         # Gene set database in GMT format
  output.directory = "",        # Directory where to store output and results (default: "")
  doc.string            = "GSEA.analysis",   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
  non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
  reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
  nperm                 = 1000,            # Number of random permutations (default: 1000)
  weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
  nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
  fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
  fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
  topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
  adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
  gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
  gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
  reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
  preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
  random.seed           = 3338,            # Random number generator seed. (default: 123456)
  perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
  fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
  replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
  save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
  OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
  use.fast.enrichment.routine = T,          # Use faster routine to compute enrichment for random permutations (default: T)
  silent = T
){
  # GSEA, we all have a love/hate relationship with it.  We love that it works
  # and is accepted.  We hate that it has so many quirks.  The most annoying
  # quirk is that we have to use this Java version because the R version is out
  # of date.  That isn't too bad.  The real problem is that we cannot specify
  # the directory where the output is going to end up.  That is torture.  If we
  # rerun an analysis the report directories can collide which makes it
  # impossible to discover the report.  To combat this I started outputting the
  # data to a UUID and then moving the directory to the correct directory when
  # done.
  
  library(uuid)
  analysis.uuid <- UUIDgenerate()
  # But GSEA has bugs.  You can't have a dash in your doc.string because it
  # chokes.  Great.
  analysis.uuid <- gsub("-", "", analysis.uuid)
  system(paste("java -cp ~/gsea/gsea2-2.2.2.jar ",
               "-Xmx4000m xtools.gsea.GseaPreranked ",
               "-gmx ", gs.db, " ",
               "-rnk ", input.rnk, " ",
               "-collapse false ",
               "-mode Max_probe ",
               "-norm meandiv ",
               "-nperm 1000 ",
               "-scoring_scheme weighted ",
               "-rpt_label ", analysis.uuid, " ",
               "-plot_top_x ", topgs, " ",
               "-rnd_seed ", random.seed, " ",
               "-save_rnd_lists false ",
               "-set_max ", gs.size.threshold.max, " ",
               "-set_min ", gs.size.threshold.min, " ", 
               "-zip_report false ",
               "-out ", output.directory, " ",
               "-gui false"), intern=silent)
  report.dir <- list.files(output.directory)
  # get the directory containing our uuid
  file.index <- max(grep(paste0("^",analysis.uuid), report.dir))
  # Switched this to hard links because rstudio doesn't follow symlinks
  final.dir <- paste(getwd(), output.directory, doc.string, sep="/")
  # I had to suppress warning for this rather than checking if it exists before
  # deleting because a broken symlink is flagged as does not exist when it does
  # in fact exist and I couldn't find a work around.
  suppressWarnings(unlink(final.dir, recursive=T))
  file.rename(paste(getwd(), output.directory, report.dir[file.index], sep="/"), 
               final.dir)
  # Return the output files
  list.files(final.dir, "gsea_report.*html", full.names=T)
}
read.gsea.results <- function(files) {
  require(XML)
  require(foreach)
  results <- foreach (file=files, .combine=rbind) %do% {
    html.in  <- htmlParse(file)
    results.table <- readHTMLTable(html.in, header=F, as.data.frame=F, which=1)
    results.table <- data.frame(lapply(results.table, type.convert, as.is=T))
    results.table <- results.table[, c(-1,-3)]
    colnames(results.table) <- c( "sig", "size", "ES", "NES", "p", "q", "fwer", "rank", "leader")
    results.table
  }
  results[order(results$ES), ]
}

find.gsea.results <- function(doc.string, output.directory="gsea") {
  report.dir <- list.files(paste(getwd(), output.directory, doc.string, sep="/"), "gsea_report.*html", full.names=T)
  report.dir
}