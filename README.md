# RNA-seq Pipeline

**Reproducible bulk RNA-seq analysis pipeline**

From raw FASTQ to differential expression, with **spike-ins**, **RUVSeq**, and **SVA/svaseq** batch correction.

Built on highly cited methods: **Salmon**, **tximport**, **DESeq2**, **RUVSeq**, **svaseq** (see [docs/CITATIONS.md](docs/CITATIONS.md)).

**Core path:**  
**QC** → Trimming → Quantification (Salmon / STAR) → Spike-in QC → **RUVSeq or SVA** → DESeq2

---

## Pipeline Overview

```
1. Quality Control
   FastQC + MultiQC + read-count warnings + spike-in recovery

2. Trimming (fastp)

3. Quantification
   Salmon (recommended) or STAR + featureCounts

4. Spike-in handling (optional)
   Recovery metrics, optional simple size factors

5. Unwanted variation / batch correction (optional)
   • RUVSeq: RUVg | RUVs | RUVr
   • SVA: svaseq (latent factors) | ComBat-seq (known batch, viz counts)
   • Diagnostics: RLE / PCA before & after

6. Differential Expression (DESeq2)
   Design can include batch, W factors, or SV factors
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

## Batch correction options

### RUVSeq
```yaml
ruv:
  enabled: true
  method: ruvg          # ruvg | ruvs | ruvr
  k: 1
  controls: spikein
```
See [docs/ruvseq.md](docs/ruvseq.md).

### SVA / svaseq
```yaml
sva:
  enabled: true
  method: svaseq        # svaseq | combat_seq
  n_sv: null            # auto-estimate
  add_to_design: true
```

```bash
Rscript scripts/sva_correction.R \
  --counts results/salmon \
  --samples data/samples.csv \
  --method svaseq \
  --out results/sva
```

See [docs/sva.md](docs/sva.md).

**Recommendation**
- Spike-ins present → RUVg  
- No controls, latent batch → svaseq  
- Known batch → put `batch` in DESeq2 design (use ComBat-seq only for visualization counts)

---

## Spike-in & QC

- Spike-ins: [docs/spikein.md](docs/spikein.md)  
- QC details: [docs/qc.md](docs/qc.md)

---

## Outputs (batch-correction related)

```
results/ruv/
  unwanted_factors.csv
  samples_with_W.csv
  RLE_*.pdf / PCA_*.pdf

results/sva/
  surrogate_variables.csv
  samples_with_SV.csv
  PCA_before.pdf / PCA_after_*.pdf
  combatseq_adjusted_counts.csv   # if method: combat_seq
```

---

## Documentation

- **[Citations (all methods)](docs/CITATIONS.md)**
- [RUVSeq (RUVg / RUVs / RUVr)](docs/ruvseq.md)
- [SVA / svaseq / ComBat-seq](docs/sva.md)
- [Spike-in controls](docs/spikein.md)
- [Quality control](docs/qc.md)
- [Design formulas & contrasts](docs/design_and_contrasts.md)

---

## Citation

If you use this pipeline, please cite the underlying tools listed in **[docs/CITATIONS.md](docs/CITATIONS.md)**. Key papers:

- Salmon — Patro et al., Nat Methods 2017  
- DESeq2 — Love et al., Genome Biol 2014  
- tximport — Soneson et al., F1000Research 2015  
- RUVSeq — Risso et al., Nat Biotechnol 2014  
- svaseq — Leek, Nucleic Acids Res 2014  

---

## License

MIT
