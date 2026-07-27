# Spike-in Controls (ERCC and custom)

This pipeline supports **spike-in controls** for quality assessment and optional normalization.

## Why use spike-ins?

- Detect global shifts in RNA content (e.g. transcriptional shut-down or amplification)
- Provide external reference points independent of endogenous genes
- Enable more accurate normalization when composition biases are strong
- QC: check whether the expected dose-response of the spike-in mix is recovered

The most widely used commercial set is the **ERCC RNA Spike-In Mix** (Thermo Fisher / Ambion).

---

## Enabling spike-ins

In `config/config.yaml`:

```yaml
spikein:
  enabled: true
  type: ERCC
  fasta: "/path/to/ERCC92.fa"
  concentration_table: "data/ERCC_Controls_Analysis.txt"   # optional
  normalization: none          # none | size_factor | ruvg
  id_prefix: "ERCC-"
```

### 1. Build an index that includes the spike-ins

**Salmon (recommended):**

```bash
# Concatenate transcriptome + ERCC
cat gencode.v44.transcripts.fa.gz <(bgzip -c ERCC92.fa) > transcriptome_plus_ERCC.fa.gz

salmon index -t transcriptome_plus_ERCC.fa.gz -i salmon_index_with_ERCC --gencode -p 8
```

**STAR:**

Include the ERCC FASTA when generating the genome index (or add them as extra sequences).

### 2. Add spike-ins to your library prep

Follow the manufacturer’s protocol (usually a fixed volume of ERCC Mix 1 or Mix 2 is added to a defined amount of total RNA).

---

## What the pipeline does with spike-ins

When `spikein.enabled: true`:

1. **QC report**
   - Fraction of reads assigned to spike-ins per sample
   - Detection rate of individual ERCC transcripts
   - Optional dose-response plot (observed vs expected abundance) if a concentration table is provided

2. **Quantification**
   - Spike-in transcripts are quantified together with endogenous transcripts (because they are in the index)

3. **Optional normalization**
   - `normalization: none` — classic DESeq2 median-of-ratios (default, usually preferred)
   - `normalization: size_factor` — estimate size factors **only from spike-in counts**
   - `normalization: ruvg` — use RUVSeq::RUVg with spike-ins as negative controls (requires `RUVSeq` package)

---

## Recommended practice

- For most bulk RNA-seq experiments **standard DESeq2 normalization is still preferred** even when spike-ins are present.
- Use spike-ins primarily for **QC and diagnosing global shifts**.
- Only switch to spike-in-based size factors or RUVg when you have strong biological reason (e.g. large changes in total RNA content per cell).

---

## ERCC resources

- ERCC92 FASTA and annotation: available from Thermo Fisher or NCBI
- Analysis file with concentrations: usually provided with the kit (`ERCC_Controls_Analysis.txt`)

Place the concentration table in `data/` and point `concentration_table` to it if you want dose-response QC plots.

---

## Output files related to spike-ins

```
results/qc/
  spikein_metrics.csv          # per-sample spike-in summary
  ercc_dose_response.pdf       # (if concentration table given)

results/deseq2/
  (standard results; size factors may be derived from spike-ins if requested)
```
