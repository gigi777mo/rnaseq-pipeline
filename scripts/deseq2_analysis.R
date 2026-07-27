#!/usr/bin/env Rscript

# Differential expression analysis with tximport + DESeq2
# Designed for Salmon quantification output

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
})

option_list <- list(
  make_option("--salmon_dir", type="character", help="Directory containing Salmon quant folders"),
  make_option("--samples", type="character", help="Sample sheet CSV"),
  make_option("--design", type="character", default="~ condition", help="DESeq2 design formula"),
  make_option("--reference", type="character", default="control", help="Reference level for condition"),
  make_option("--fdr", type="double", default=0.05),
  make_option("--lfc", type="double", default=1.0),
  make_option("--out", type="character", default="results/deseq2")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$out, recursive=TRUE, showWarnings=FALSE)

# ---------- Load sample sheet ----------
samples <- read.csv(opt$samples, stringsAsFactors=FALSE)
rownames(samples) <- samples$sample

# Build file list for tximport
files <- file.path(opt$salmon_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

missing <- !file.exists(files)
if (any(missing)) {
  stop("Missing quant.sf files for: ", paste(names(files)[missing], collapse=", "))
}

# ---------- tximport ----------
# Note: for gene-level analysis we need a tx2gene map.
# If you have a GTF, create tx2gene beforehand or use tximeta.
# Here we assume the user provides transcript-level or has pre-aggregated.
# For simplicity in this starter pipeline we import and let the user
# supply a proper tx2gene if needed. Many pipelines use tximeta for this.

message("Importing Salmon quantifications...")
txi <- tximport(files, type="salmon", txOut=TRUE)   # transcript-level first

# Simple gene-level aggregation if transcript names contain gene info (GENCODE style)
# GENCODE: ENST...|ENSG...|...
# For production use, generate a proper tx2gene table from GTF.

message("Running DESeq2...")

# Make sure condition is a factor with the desired reference level
samples$condition <- factor(samples$condition)
if (opt$reference %in% levels(samples$condition)) {
  samples$condition <- relevel(samples$condition, ref=opt$reference)
}

dds <- DESeqDataSetFromTximport(txi, colData=samples, design=as.formula(opt$design))

# Pre-filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# Results
res <- results(dds, alpha=opt$fdr)
res <- lfcShrink(dds, coef=resultsNames(dds)[length(resultsNames(dds))], type="apeglm")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, file.path(opt$out, "results.csv"), row.names=FALSE)

sig <- subset(res_df, padj < opt$fdr & abs(log2FoldChange) >= opt$lfc)
write.csv(sig, file.path(opt$out, "significant.csv"), row.names=FALSE)

# Normalized counts
norm <- counts(dds, normalized=TRUE)
write.csv(norm, file.path(opt$out, "normalized_counts.csv"))

# ---------- Plots ----------
# PCA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA")
ggsave(file.path(opt$out, "pca.pdf"), p, width=6, height=5)

# Volcano
pdf(file.path(opt$out, "volcano.pdf"), width=8, height=7)
print(EnhancedVolcano(res_df,
                      lab = res_df$gene,
                      x = "log2FoldChange",
                      y = "padj",
                      pCutoff = opt$fdr,
                      FCcutoff = opt$lfc,
                      title = "Volcano plot"))
dev.off()

message("DESeq2 analysis complete. Results written to ", opt$out)
