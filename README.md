# RNA-seq Pipeline

**Reproducible bulk RNA-seq analysis pipeline**

From raw FASTQ files to differential expression results.

**Core path:**  
**QC** (FastQC / MultiQC + read-count checks) → Trimming (fastp) → Quantification (**Salmon** or STAR) → **Spike-in QC & optional normalization** → Differential expression (**DESeq2**)

---

## Pipeline Overview

```
1. Quality Control (enhanced)
   • FastQC on raw (and optionally trimmed) reads
   • MultiQC HTML report
   • Per-sample read-count warnings
   • Spike-in recovery metrics (when enabled)

2. Adapter & quality trimming
   fastp

3. Quantification
   A. Salmon (alignment-free, recommended)
   B. STAR + featureCounts

4. Spike-in handling (optional)
   • ERCC or custom spike-ins
   • Recovery & dose-response QC
   • Optional size-factor or RUVg normalization

5. Differential Expression
   tximport → DESeq2 (with LFC shrinkage)
   PCA, volcano, results tables
```

---

## Quick Start

```bash
git clone https://github.com/gigi777mo/rnaseq-pipeline.git
cd rnaseq-pipeline

conda env create -f environment.yml
conda activate rnaseq

# Edit config/config.yaml and data/samples.csv
# Place paired-end FASTQs in data/raw/

snakemake --cores 8 --use-conda
```

Or run steps manually:

```bash
python scripts/run_qc.py --input data/raw --out results/qc --stage raw
python scripts/run_trim.py --input data/raw --out results/trimmed
python scripts/run_salmon.py --input results/trimmed --index /path/to/index --out results/salmon
python scripts/run_qc.py --input results/trimmed --out results/qc --stage trimmed \
    --salmon-dir results/salmon --samples data/samples.csv   # includes spike-in QC if present
Rscript scripts/deseq2_analysis.R --salmon_dir results/salmon --samples data/samples.csv --out results/deseq2
```

---

## Quality Control (enhanced)

The QC module now provides:

| Check | Description |
|-------|-------------|
| FastQC | Per-base quality, adapter content, duplication, GC, etc. |
| MultiQC | Single interactive HTML report aggregating all samples |
| Read-count warning | Flags samples below a configurable minimum number of reads |
| Spike-in metrics | Fraction of reads on spike-ins, number detected, recovery plot |

Configuration (`config/config.yaml`):

```yaml
qc:
  run_fastqc: true
  run_multiqc: true
  min_reads: 1000000
  min_mapping_rate: 0.5
```

---

## Spike-in Support (ERCC or custom)

Enable in `config/config.yaml`:

```yaml
spikein:
  enabled: true
  type: ERCC
  fasta: "/path/to/ERCC92.fa"
  concentration_table: "data/ERCC_Controls_Analysis.txt"   # optional
  normalization: none          # none | size_factor | ruvg
  id_prefix: "ERCC-"
```

**What you get:**

- Spike-in transcripts quantified together with endogenous genes (index must contain them)
- `results/qc/spikein_metrics.csv` – per-sample recovery summary
- `results/qc/spikein_fraction.pdf` – visual QC
- Optional dose-response plot when a concentration table is supplied
- Optional alternative normalization (spike-in size factors or RUVg)

Full documentation: **[docs/spikein.md](docs/spikein.md)**

> **Recommendation:** Use spike-ins primarily for QC. Standard DESeq2 median-of-ratios normalization is still preferred for most experiments unless you expect large global changes in RNA content.

---

## Configuration highlights

```yaml
genome: human
quantifier: salmon              # or star

qc:
  min_reads: 1000000

spikein:
  enabled: false                # turn on when you have ERCC/custom spike-ins
  normalization: none

deseq2:
  design: "~ condition"
  reference_level: "control"
  fdr_cutoff: 0.05
  lfc_cutoff: 1.0
```

---

## Outputs

```
results/
├── qc/
│   ├── multiqc_report.html
│   ├── spikein_metrics.csv      # when spike-ins enabled
│   └── spikein_fraction.pdf
├── trimmed/
├── salmon/   (or star/)
├── deseq2/
│   ├── results.csv
│   ├── significant.csv
│   ├── pca.pdf
│   ├── volcano.pdf
│   └── normalized_counts.csv
└── logs/
```

---

## Building an index with spike-ins

```bash
# Example for Salmon + ERCC
cat gencode.v44.transcripts.fa.gz <(bgzip -c ERCC92.fa) > tx_plus_ERCC.fa.gz
salmon index -t tx_plus_ERCC.fa.gz -i salmon_index_ERCC --gencode -p 8
```

Then set `salmon_index` in the config to this new index.

---

## Tips

- Always inspect the MultiQC report before differential expression.
- ≥3 biological replicates per group is strongly recommended.
- Include `batch` in the design formula when appropriate: `~ batch + condition`.
- Spike-ins are excellent for diagnosing global effects; use them for normalization only when biologically justified.

---

## Citation

- Salmon — Patro et al., Nat Methods 2017
- DESeq2 — Love et al., Genome Biol 2014
- STAR — Dobin et al., Bioinformatics 2013
- fastp — Chen et al., Bioinformatics 2018
- tximport — Soneson et al., F1000Research 2015
- ERCC spike-ins — External RNA Controls Consortium

---

## License

MIT
