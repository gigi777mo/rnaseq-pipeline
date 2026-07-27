#!/usr/bin/env Rscript

# RUVSeq normalization / factor estimation
# Supports RUVg, RUVs, and RUVr with diagnostic plots

suppressPackageStartupMessages({
  library(optparse)
  library(RUVSeq)
  library(DESeq2)
  library(ggplot2)
  library(matrixStats)
})

option_list <- list(
  make_option("--counts", type="character",
              help="Count matrix CSV (genes x samples) or Salmon directory"),
  make_option("--samples", type="character", help="Sample sheet CSV"),
  make_option("--method", type="character", default="ruvg",
              help="ruvg | ruvs | ruvr"),
  make_option("--k", type="integer", default=1),
  make_option("--controls", type="character", default="spikein",
              help="spikein | empirical | custom"),
  make_option("--spikein-prefix", type="character", default="ERCC-"),
  make_option("--custom-controls", type="character", default="",
              help="File with control gene IDs (one per line)"),
  make_option("--replicate-column", type="character", default="",
              help="Column in sample sheet for RUVs groups"),
  make_option("--design", type="character", default="~ condition"),
  make_option("--out", type="character", default="results/ruv")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out, recursive=TRUE, showWarnings=FALSE)

# ---------- Load data ----------
samples <- read.csv(opt$samples, stringsAsFactors=FALSE)
rownames(samples) <- samples$sample

if (dir.exists(opt$counts)) {
  # Assume Salmon quant folders – build a simple count matrix from NumReads
  # (For production use prefer tximport → gene-level counts)
  message("Reading Salmon quant.sf files...")
  files <- file.path(opt$counts, samples$sample, "quant.sf")
  names(files) <- samples$sample
  mats <- lapply(files, function(f) {
    q <- read.delim(f, stringsAsFactors=FALSE)
    setNames(q$NumReads, q$Name)
  })
  genes <- unique(unlist(lapply(mats, names)))
  counts <- matrix(0, nrow=length(genes), ncol=length(mats),
                   dimnames=list(genes, names(mats)))
  for (s in names(mats)) counts[names(mats[[s]]), s] <- mats[[s]]
} else {
  counts <- as.matrix(read.csv(opt$counts, row.names=1, check.names=FALSE))
  counts <- counts[, samples$sample, drop=FALSE]
}

# Filter very low counts
keep <- rowSums(counts) >= 10
counts <- counts[keep, ]
message("Genes/transcripts retained: ", nrow(counts))

# Create SeqExpressionSet
set <- newSeqExpressionSet(counts, phenoData=data.frame(samples, row.names=samples$sample))

# ---------- Helper: diagnostics ----------
save_diagnostics <- function(set_obj, tag) {
  # RLE
  pdf(file.path(opt$out, paste0("RLE_", tag, ".pdf")), width=8, height=5)
  plotRLE(set_obj, outline=FALSE, ylim=c(-4, 4), main=paste("RLE -", tag))
  dev.off()

  # PCA
  pdf(file.path(opt$out, paste0("PCA_", tag, ".pdf")), width=6, height=5)
  plotPCA(set_obj, main=paste("PCA -", tag), col=as.numeric(as.factor(pData(set_obj)$condition)))
  dev.off()
}

save_diagnostics(set, "before_RUV")

# ---------- Run chosen method ----------
method <- tolower(opt$method)
k <- opt$k

if (method == "ruvg") {
  message("Running RUVg (k=", k, ")...")

  if (opt$controls == "spikein") {
    cIdx <- grepl(opt$`spikein-prefix`, rownames(counts))
    if (sum(cIdx) < 5) stop("Too few spike-in controls found with prefix: ", opt$`spikein-prefix`)
    message("Using ", sum(cIdx), " spike-in controls")
  } else if (opt$controls == "custom") {
    ids <- readLines(opt$`custom-controls`)
    cIdx <- rownames(counts) %in% ids
    message("Using ", sum(cIdx), " custom controls")
  } else if (opt$controls == "empirical") {
    # First-pass edgeR-like ranking via simple variance / mean filter as placeholder
    # In practice a proper DE test is better; here we take lowest-variance genes
    message("Selecting empirical controls (lowest variance genes)...")
    vars <- rowVars(log1p(counts))
    n_ctrl <- min(5000, floor(nrow(counts) * 0.5))
    cIdx <- rank(vars) <= n_ctrl
  } else {
    stop("Unknown controls mode: ", opt$controls)
  }

  set_ruv <- RUVg(set, cIdx, k=k)

} else if (method == "ruvs") {
  message("Running RUVs (k=", k, ")...")
  if (!nzchar(opt$`replicate-column`) || !opt$`replicate-column` %in% colnames(samples)) {
    stop("RUVs requires --replicate-column pointing to a column in the sample sheet")
  }
  groups <- matrix(samples[[opt$`replicate-column`]], ncol=1)
  rownames(groups) <- samples$sample
  # RUVs expects a makeGroups-style list; simplify via differences within groups
  set_ruv <- RUVs(set, cIdx=rownames(counts), k=k, scIdx=makeGroups(samples[[opt$`replicate-column`]]))

} else if (method == "ruvr") {
  message("Running RUVr (k=", k, ")...")
  # Design matrix for first-pass residuals
  design_formula <- as.formula(opt$design)
  mm <- model.matrix(design_formula, data=samples)
  # Use deviance residuals via edgeR-style approach inside RUVr
  set_ruv <- RUVr(set, cIdx=rownames(counts), k=k, residuals=residuals(glm.fit(...))) # fallback below

  # Practical implementation: compute deviance residuals with a simple NB GLM approximation
  # For robustness we use the built-in approach of passing upper-quartile normalized residuals
  library(edgeR)
  dge <- DGEList(counts=counts)
  dge <- calcNormFactors(dge, method="upperquartile")
  design <- model.matrix(as.formula(opt$design), data=samples)
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  res <- residuals(fit, type="deviance")
  set_ruv <- RUVr(set, cIdx=rownames(counts), k=k, residuals=res)

} else {
  stop("Unknown method: ", method)
}

save_diagnostics(set_ruv, paste0("after_", method))

# ---------- Save outputs ----------
# Estimated factors
w_factors <- pData(set_ruv)[, grep("^W_", colnames(pData(set_ruv))), drop=FALSE]
write.csv(w_factors, file.path(opt$out, "unwanted_factors.csv"))

# Normalized counts
norm <- normCounts(set_ruv)
write.csv(norm, file.path(opt$out, "normalized_counts.csv"))

# Updated sample sheet with W factors for DESeq2
samples_out <- cbind(samples, w_factors)
write.csv(samples_out, file.path(opt$out, "samples_with_W.csv"), row.names=FALSE)

message("RUV finished. Outputs written to ", opt$out)
message("Add the W_* columns to your DESeq2 design, e.g.: ~ W_1 + condition")
