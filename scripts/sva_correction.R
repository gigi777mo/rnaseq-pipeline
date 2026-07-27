#!/usr/bin/env Rscript

# SVA / svaseq / ComBat-seq batch correction
# Estimates surrogate variables or applies ComBat-seq and writes diagnostics

suppressPackageStartupMessages({
  library(optparse)
  library(sva)
  library(DESeq2)
  library(ggplot2)
})

option_list <- list(
  make_option("--counts", type="character",
              help="Count matrix CSV (genes x samples) or Salmon quant directory"),
  make_option("--samples", type="character", help="Sample sheet CSV"),
  make_option("--method", type="character", default="svaseq",
              help="svaseq | combat_seq"),
  make_option("--design", type="character", default="~ condition",
              help="Full model formula (biological factors)"),
  make_option("--n-sv", type="integer", default=NA,
              help="Number of SVs (NA = estimate with num.sv)"),
  make_option("--num-sv-method", type="character", default="be",
              help="be | leek"),
  make_option("--batch-column", type="character", default="batch",
              help="Column name for known batch (ComBat-seq)"),
  make_option("--out", type="character", default="results/sva")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out, recursive=TRUE, showWarnings=FALSE)

# ---------- Load samples & counts ----------
samples <- read.csv(opt$samples, stringsAsFactors=FALSE)
rownames(samples) <- samples$sample

if (dir.exists(opt$counts)) {
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

keep <- rowSums(counts) >= 10
counts <- counts[keep, ]
message("Features retained: ", nrow(counts))

# Ensure condition factor exists
if (!"condition" %in% colnames(samples)) {
  stop("Sample sheet must contain a 'condition' column")
}
samples$condition <- factor(samples$condition)

# ---------- Diagnostics helper ----------
plot_pca <- function(mat, coldata, title, file) {
  # Simple log-CPM style PCA for visualization
  logmat <- log1p(mat)
  pc <- prcomp(t(logmat), center=TRUE, scale.=FALSE)
  percent <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], condition=coldata$condition)
  if ("batch" %in% colnames(coldata)) df$batch <- coldata$batch
  p <- ggplot(df, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percent[1], "%")) +
    ylab(paste0("PC2: ", percent[2], "%")) +
    theme_bw() + ggtitle(title)
  ggsave(file, p, width=6, height=5)
}

plot_pca(counts, samples, "PCA before correction", file.path(opt$out, "PCA_before.pdf"))

method <- tolower(opt$method)

if (method == "svaseq") {
  message("Running svaseq...")
  # Moderated log is handled inside svaseq; pass counts
  mod  <- model.matrix(as.formula(opt$design), data=samples)
  mod0 <- model.matrix(~ 1, data=samples)

  n.sv <- opt$`n-sv`
  if (is.na(n.sv)) {
    n.sv <- num.sv(counts, mod, method=opt$`num-sv-method`)
    message("Estimated number of SVs: ", n.sv)
  }
  if (n.sv < 1) {
    message("No surrogate variables estimated (n.sv = 0). Writing empty result.")
    write.csv(data.frame(row.names=samples$sample), file.path(opt$out, "surrogate_variables.csv"))
    write.csv(samples, file.path(opt$out, "samples_with_SV.csv"), row.names=FALSE)
    quit(save="no", status=0)
  }

  svobj <- svaseq(counts, mod, mod0, n.sv=n.sv)
  sv <- as.data.frame(svobj$sv)
  colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
  rownames(sv) <- samples$sample

  write.csv(sv, file.path(opt$out, "surrogate_variables.csv"))
  samples_out <- cbind(samples, sv)
  write.csv(samples_out, file.path(opt$out, "samples_with_SV.csv"), row.names=FALSE)

  # Residual-style PCA after regressing SVs (approximate visualization)
  # For true DE, user should add SVs to design; here we show residualization for QC only
  design_sv <- model.matrix(~ . , data=sv)
  resid <- t(lm(t(log1p(counts)) ~ design_sv)$residuals)
  plot_pca(exp(resid) - 1, samples, "PCA after SV residualization (viz only)",
           file.path(opt$out, "PCA_after_svaseq.pdf"))

  message("svaseq complete. Add SV columns to DESeq2 design, e.g. ~ SV1 + condition")

} else if (method == "combat_seq") {
  message("Running ComBat-seq...")
  if (!opt$`batch-column` %in% colnames(samples)) {
    stop("ComBat-seq requires a batch column: ", opt$`batch-column`)
  }
  batch <- samples[[opt$`batch-column`]]
  group <- samples$condition

  adjusted <- ComBat_seq(counts, batch=batch, group=group, full_mod=TRUE)
  write.csv(adjusted, file.path(opt$out, "combatseq_adjusted_counts.csv"))

  plot_pca(adjusted, samples, "PCA after ComBat-seq (use for viz only)",
           file.path(opt$out, "PCA_after_combatseq.pdf"))

  message("ComBat-seq complete.")
  message("WARNING: Use adjusted counts for visualization only.")
  message("For DESeq2, keep original counts and put batch in the design formula.")

} else {
  stop("Unknown method: ", method, " (use svaseq or combat_seq)")
}

message("Outputs written to ", opt$out)
