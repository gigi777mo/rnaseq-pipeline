# RUVSeq Normalization Methods

This pipeline supports the three main procedures from the Bioconductor package **RUVSeq** (Risso et al., Nat Biotechnol 2014) for removing unwanted variation from RNA-seq counts.

---

## Overview of the three methods

| Method | Uses | Best when | Main risk |
|--------|------|-----------|-----------|
| **RUVg** | Negative control **genes** (spike-ins, housekeepers, empirical) | You have ERCC or trusted controls | Poor controls → over/under-correction |
| **RUVs** | Negative control **samples** / technical replicates | You have proper replicates | Needs replicate design |
| **RUVr** | **Residuals** from a first-pass GLM | No controls available | Can remove biology if model is wrong |

All three estimate factors of unwanted variation (*W*) that can be:
- regressed out to produce normalized counts, or
- added as covariates in the DESeq2 / edgeR design (`~ W_1 + condition`).

---

## 1. RUVg — control genes

**Assumption:** A set of genes is not differentially expressed with respect to the biological factor of interest, yet still captures the unwanted variation.

**Typical controls:**
- ERCC spike-ins (recommended when available)
- Housekeeping genes
- Empirical controls (least DE genes from a first-pass analysis)

**Pipeline usage:**
```yaml
spikein:
  enabled: true
  normalization: ruvg
  ruvg_k: 1
  id_prefix: "ERCC-"
```

Or more generally in the RUV script:
```r
RUVg(counts, cIdx = ercc_ids, k = 1)
```

**Notes:**
- Start with `k = 1`. Increase only while RLE and PCA diagnostics improve.
- Spike-ins are convenient but not always perfect; check that they behave consistently across samples.

---

## 2. RUVs — control samples / replicates

**Assumption:** You have technical replicates (or negative-control samples) for which the biological covariates are constant.

**How it works:** Factor analysis is performed on centered counts within replicate groups.

**Pipeline usage:**
```yaml
ruv:
  method: ruvs
  k: 1
  # sample sheet must contain a column that defines replicate groups
  # e.g. patient_id or library_id
```

RUVs is more robust to the exact choice of control genes than RUVg and can even use all genes as controls in many cases.

---

## 3. RUVr — residuals

**Assumption:** None about controls. Uses residuals from a first-pass GLM of the counts on the known design factors.

**How it works:**
1. Fit GLM ~ design (e.g. condition)
2. Compute residuals (usually deviance residuals)
3. Factor analysis on residuals → estimate *W*
4. Include *W* in a second-pass model or produce adjusted counts

**Pipeline usage:**
```yaml
ruv:
  method: ruvr
  k: 1
```

Useful when spike-ins and replicates are unavailable. Be cautious if biology is correlated with technical batches.

---

## Recommended workflow in this pipeline

1. Run standard QC + quantification (Salmon/STAR).
2. If spike-ins are present → try **RUVg** first (`k = 1`).
3. Inspect diagnostics:
   - Relative Log Expression (RLE) plots before/after
   - PCA before/after (unwanted variation should leave the leading PCs)
   - p-value histogram of the DE test (should be flat under the null + enrichment of small p-values)
4. If no spike-ins but technical replicates exist → **RUVs**.
5. Otherwise → **RUVr** or empirical-control RUVg.
6. Prefer adding the estimated *W* factors as covariates in DESeq2 rather than only using the normalized counts, when possible:
   ```r
   design(dds) <- ~ W_1 + condition
   ```

---

## Choosing *k* (number of factors)

- Start with `k = 1`.
- Increase gradually (2, 3 …) while watching RLE and PCA.
- Stop when additional factors no longer improve sample clustering by the biological variable of interest or begin to mix biology into *W*.
- Over-correcting (too large *k*) can remove true biological signal.

---

## Relation to simple spike-in size factors

Estimating size factors **only** from ERCC counts (`normalization: size_factor`) is a simpler alternative. The original RUVSeq paper showed that plain ERCC global scaling is often less effective than RUVg. Prefer RUVg when using spike-ins for normalization.

---

## Key references

- Risso D, Ngai J, Speed TP, Dudoit S. *Normalization of RNA-seq data using factor analysis of control genes or samples.* Nat Biotechnol. 2014.
- Gagnon-Bartsch JA, Speed TP. *Using control genes to correct for unwanted variation in microarray data.* Biostatistics. 2012.
- Bioconductor package: https://bioconductor.org/packages/RUVSeq
