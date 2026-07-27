# Quality Control in this pipeline

## Stages

1. **Raw-read QC** (FastQC)
   - Per-base sequence quality
   - Adapter content
   - Per-sequence GC content
   - Sequence duplication levels
   - Overrepresented sequences

2. **Trimmed-read QC** (optional second FastQC pass)
   - Confirms adapters were removed and quality improved

3. **MultiQC**
   - Aggregates all FastQC (and other tool) reports into one interactive HTML page

4. **Read-count checks**
   - Warns when a sample falls below `qc.min_reads` (default 1 million)

5. **Spike-in QC** (when enabled)
   - Fraction of total reads assigned to spike-ins
   - Number of spike-in species detected
   - Per-sample recovery barplot
   - Optional observed-vs-expected dose-response (ERCC)

## Recommended thresholds (starting points)

| Metric | Typical warning level |
|--------|-----------------------|
| Total reads (PE) | < 10–20 million for bulk mRNA-seq (experiment-dependent) |
| Config `min_reads` | 1 000 000 (very low; raise for production) |
| Spike-in fraction | Usually 0.5–5 % depending on how much was added |
| Mapping / quantification rate | > 70–80 % for good libraries |

Adjust `qc.min_reads` and `qc.min_mapping_rate` in `config/config.yaml` to match your experimental design.

## Files produced

```
results/qc/
  raw/ or trimmed/          # FastQC html/zip
  multiqc_report.html
  spikein_metrics.csv       # if spike-ins used
  spikein_fraction.pdf
```
