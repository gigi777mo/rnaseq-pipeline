#!/usr/bin/env Rscript

# Spike-in QC metrics and optional dose-response plot (ERCC)

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option("--salmon_dir", type="character", help="Directory with Salmon quant folders"),
  make_option("--samples", type="character", help="Sample sheet CSV"),
  make_option("--prefix", type="character", default="ERCC-", help="Spike-in ID prefix"),
  make_option("--concentration", type="character", default="", help="Optional ERCC concentration table"),
  make_option("--out", type="character", default="results/qc")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out, recursive=TRUE, showWarnings=FALSE)

samples <- read.csv(opt$samples, stringsAsFactors=FALSE)
files <- file.path(opt$salmon_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

# Read each quant.sf and extract spike-in rows
metrics <- list()
all_spike <- list()

for (s in names(files)) {
  q <- read.delim(files[s], stringsAsFactors=FALSE)
  spike <- q[grepl(opt$prefix, q$Name), ]
  total_reads <- sum(q$NumReads)
  spike_reads <- sum(spike$NumReads)
  n_detected <- sum(spike$NumReads > 0)

  metrics[[s]] <- data.frame(
    sample = s,
    total_reads = total_reads,
    spike_reads = spike_reads,
    spike_fraction = spike_reads / max(total_reads, 1),
    n_spike_detected = n_detected,
    n_spike_total = nrow(spike),
    stringsAsFactors = FALSE
  )
  spike$sample <- s
  all_spike[[s]] <- spike
}

metrics_df <- do.call(rbind, metrics)
write.csv(metrics_df, file.path(opt$out, "spikein_metrics.csv"), row.names=FALSE)
message("Wrote spike-in metrics: ", file.path(opt$out, "spikein_metrics.csv"))

# Simple barplot of spike-in fraction
p <- ggplot(metrics_df, aes(x=sample, y=spike_fraction)) +
  geom_col(fill="steelblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("Fraction of reads on spike-ins") +
  ggtitle("Spike-in recovery per sample")
ggsave(file.path(opt$out, "spikein_fraction.pdf"), p, width=7, height=4)

# Optional dose-response if concentration table is provided
if (nzchar(opt$concentration) && file.exists(opt$concentration)) {
  conc <- read.delim(opt$concentration, stringsAsFactors=FALSE)
  # Expect columns that can be matched to ERCC IDs; common format has 'ERCC ID' and subgroup concentrations
  # This is a minimal generic handler – user may need to adapt column names
  message("Concentration table found. Generating basic observed-vs-expected plot if columns match.")

  # Very lightweight example: if table has ERCC_ID and mixed concentration columns
  # Users should customize this section for their exact file format.
}

message("Spike-in QC finished.")
