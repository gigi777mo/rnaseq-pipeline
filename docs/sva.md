# SVA / svaseq Batch Correction

This pipeline supports **Surrogate Variable Analysis** via the Bioconductor `sva` package, complementing the existing RUVSeq methods.

---

## When to use SVA

| Situation | Recommended tool |
|-----------|------------------|
| Batch **known**, not confounded | Put `batch` in DESeq2 design |
| Batch known, need corrected counts for plots | ComBat-seq (visualization only) |
| Batch **unknown**, have spike-ins / controls | **RUVg** (preferred) or supervised SVA |
| Batch unknown, have replicates | RUVs or **svaseq** |
| Batch unknown, **no controls** | **svaseq** or RUVr |

**Rule of thumb:** Prefer adding surrogate variables (or RUV *W* factors) to the DESeq2 design formula. Do **not** feed ComBat-seq–adjusted counts into DESeq2 together with a batch term (double-correction).

---

## svaseq (RNA-seq specific)

`svaseq` is the sequencing-adapted version of SVA (Leek 2014):

1. Applies a moderated log transform to counts/FPKM.
2. Iteratively identifies empirical control genes.
3. Estimates surrogate variables (SVs) from those genes.
4. SVs are added as covariates: `~ SV1 + SV2 + condition`.

**Key parameters:**
- `n.sv` — number of surrogate variables (estimated by `num.sv` or set manually)
- `method` — `"irw"` (default), `"two-step"`, or `"supervised"`
- Supervised mode can use known control probes (similar spirit to RUVg)

---

## ComBat-seq

For **known** batches when you need an adjusted count matrix (e.g. for a collaborator or non-model-based tools):

```r
adjusted <- ComBat_seq(counts, batch = batch, group = condition, full_mod = TRUE)
```

Use adjusted counts for PCA/heatmaps only. Keep original counts for DESeq2.

---

## Pipeline usage

```yaml
sva:
  enabled: true
  method: svaseq          # svaseq | combat_seq
  n_sv: null              # null = estimate with num.sv; or set an integer
  num_sv_method: be       # be | leek
  add_to_design: true
  diagnostics: true
```

Standalone script:

```bash
Rscript scripts/sva_correction.R \
  --counts results/salmon \
  --samples data/samples.csv \
  --method svaseq \
  --out results/sva
```

**Outputs:**
- `surrogate_variables.csv`
- `samples_with_SV.csv` (ready for DESeq2)
- PCA before/after (when diagnostics enabled)

---

## SVA vs RUVSeq (short)

- **RUVg** — best when you have good negative controls (ERCC, housekeepers).
- **svaseq** — strong default when controls are absent; estimates latent factors from the data.
- **RUVr** — residual-based alternative to svaseq.
- **ComBat-seq** — known batch, adjusted counts for visualization.

Both SVA and RUV factors should be included in the DE design rather than used only to rewrite the count matrix for differential expression.

---

## References

- Leek JT. *svaseq: removing batch effects and other unwanted noise from sequencing data.* Nucleic Acids Res. 2014.
- Leek JT, Storey JD. *Capturing heterogeneity in gene expression studies by surrogate variable analysis.* PLoS Genet. 2007.
- Johnson WE et al. *Adjusting batch effects in microarray expression data using empirical Bayes methods.* Biostatistics. 2007. (ComBat)
- Zhang Y et al. *ComBat-seq* (negative-binomial extension for RNA-seq counts).
