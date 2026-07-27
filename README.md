# RNA-seq Pipeline

**Reproducible bulk RNA-seq analysis pipeline**

From raw FASTQ files to differential expression results, with full support for **spike-ins** and **RUVSeq** normalization.

**Core path:**  
**QC** → Trimming → Quantification (Salmon / STAR) → **Spike-in QC** → **RUVSeq (RUVg / RUVs / RUVr)** → DESeq2

---

## Pipeline Overview

```
1. Quality Control
   FastQC + MultiQC + read-count warnings + spike-in recovery metrics

2. Trimming (fastp)

3. Quantification
   Salmon (recommended) or STAR + featureCounts

4. Spike-in handling (optional)
   Recovery metrics, dose-response, simple size-factor option

5. Remove Unwanted Variation (optional)
   RUVg (control genes / spike-ins)
   RUVs (replicate samples)
   RUVr (residuals)
   + RLE & PCA diagnostics before/after

6. Differential Expression (DESeq2)
   Can include estimated W factors in the design
```

---

## Quick Start

```bash
git clone https://github.com/gigi777mo/rnaseq-pipeline.git
cd rnaseq-pipeline
conda env create -f environment.yml
conda activate rnaseq

# Edit config/config.yaml and data/samples.csv
# Place FASTQs in data/raw/

snakemake --cores 8 --use-conda
```

---

## RUVSeq Support (RUVg / RUVs / RUVr)

Enable in `config/config.yaml`:

```yaml
ruv:
  enabled: true
  method: ruvg          # ruvg | ruvs | ruvr
  k: 1
  controls: spikein     # spikein | empirical | custom
  add_to_design: true
  diagnostics: true
```

Or run the standalone script:

```bash
Rscript scripts/ruv_normalization.R \
  --counts results/salmon \
  --samples data/samples.csv \
  --method ruvg \
  --k 1 \
  --controls spikein \
  --spikein-prefix ERCC- \
  --out results/ruv
```

**Outputs:**
- `unwanted_factors.csv` – estimated *W* factors
- `normalized_counts.csv`
- `samples_with_W.csv` – ready to use in DESeq2 (`~ W_1 + condition`)
- RLE and PCA plots before/after correction

Full explanation of the three methods: **[docs/ruvseq.md](docs/ruvseq.md)**

---

## Spike-in Support

```yaml
spikein:
  enabled: true
  type: ERCC
  id_prefix: "ERCC-"
  normalization: none     # none | size_factor (simple alternative to RUVg)
```

See **[docs/spikein.md](docs/spikein.md)**.

---

## Quality Control

Enhanced QC includes FastQC, MultiQC, read-count warnings, and spike-in recovery metrics.  
Details: **[docs/qc.md](docs/qc.md)**

---

## Configuration (key sections)

```yaml
qc:
  min_reads: 1000000

spikein:
  enabled: false
  normalization: none

ruv:
  enabled: false
  method: ruvg              # ruvg | ruvs | ruvr
  k: 1
  controls: spikein
  add_to_design: true

deseq2:
  design: "~ condition"     # becomes "~ W_1 + condition" when RUV adds factors
  fdr_cutoff: 0.05
  lfc_cutoff: 1.0
```

---

## Outputs

```
results/
├── qc/
│   ├── multiqc_report.html
│   ├── spikein_metrics.csv
│   └── spikein_fraction.pdf
├── ruv/                    # when RUV enabled
│   ├── unwanted_factors.csv
│   ├── samples_with_W.csv
│   ├── normalized_counts.csv
│   ├── RLE_before_RUV.pdf / RLE_after_*.pdf
│   └── PCA_before_RUV.pdf / PCA_after_*.pdf
├── deseq2/
│   ├── results.csv
│   ├── significant.csv
│   ├── pca.pdf / volcano.pdf
│   └── normalized_counts.csv
└── ...
```

---

## Tips

- Always inspect MultiQC + RLE/PCA before trusting DE results.
- With spike-ins, prefer **RUVg** over simple spike-in size factors.
- Start with `k = 1` and increase only while diagnostics improve.
- Adding *W* factors to the DESeq2 design is usually better than only using RUV-normalized counts.
- ≥3 biological replicates per group recommended.

---

## Documentation

- [RUVSeq methods (RUVg / RUVs / RUVr)](docs/ruvseq.md)
- [Spike-in controls](docs/spikein.md)
- [Quality control](docs/qc.md)
- [Design formulas & contrasts](docs/design_and_contrasts.md)

---

## Citation

- Salmon — Patro et al., Nat Methods 2017
- DESeq2 — Love et al., Genome Biol 2014
- RUVSeq — Risso et al., Nat Biotechnol 2014
- STAR, fastp, tximport, ERCC spike-ins — see respective papers

---

## License

MIT
