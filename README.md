# RNA-seq Pipeline

**Reproducible bulk RNA-seq analysis pipeline**

From raw FASTQ files to differential expression results.

**Core path:**  
QC (FastQC/MultiQC) → Trimming (fastp) → Quantification (**Salmon** recommended or STAR + featureCounts) → Differential expression (**DESeq2**)

Designed to be modular, well-documented, and easy to adapt for human or mouse experiments.

---

## Pipeline Overview

```
1. Quality Control
   FastQC → MultiQC report

2. Adapter & quality trimming
   fastp (fast and accurate)

3. Quantification (choose one)
   A. Salmon (alignment-free, recommended – fast + isoform-aware)
   B. STAR alignment + featureCounts (classic alignment-based)

4. Differential Expression
   tximport (if Salmon) → DESeq2
   Results tables, MA/volcano plots, PCA

5. Optional downstream
   Pathway enrichment, custom contrasts, batch correction
```

---

## Quick Start

### 1. Install

```bash
git clone https://github.com/gigi777mo/rnaseq-pipeline.git
cd rnaseq-pipeline

# Create environment
conda env create -f environment.yml
conda activate rnaseq

# Optional: install Snakemake if you want the full workflow manager
conda install -c bioconda -c conda-forge snakemake
```

### 2. Prepare your data

```
data/
├── raw/                  # place paired-end FASTQ files here
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── ...
├── samples.csv           # sample metadata (see example)
└── config/
    └── config.yaml         # pipeline parameters
```

Edit `data/samples.csv` and `config/config.yaml`.

### 3. Run

**Option A – Snakemake (recommended for full automation)**
```bash
snakemake --cores 8 --use-conda
```

**Option B – Step-by-step scripts**
```bash
# QC
python scripts/run_qc.py --input data/raw --out results/qc

# Trim
python scripts/run_trim.py --input data/raw --out results/trimmed

# Quantify with Salmon
python scripts/run_salmon.py --input results/trimmed --index /path/to/salmon_index --out results/salmon

# Differential expression
Rscript scripts/deseq2_analysis.R --counts results/salmon --samples data/samples.csv --out results/deseq2
```

---

## Configuration

### `data/samples.csv`

```csv
sample,condition,batch,replicate
control_1,control,1,1
control_2,control,1,2
treated_1,treated,1,1
treated_2,treated,1,2
```

### `config/config.yaml` (key options)

```yaml
genome: human                 # human or mouse
quantifier: salmon            # salmon or star
fdr_cutoff: 0.05
lfc_cutoff: 1.0               # log2 fold-change threshold for plots
design: "~ condition"         # DESeq2 design formula

resources:
  threads: 8
  memory_gb: 32
```

---

## Quantification Choices

| Method | Pros | Cons | When to use |
|--------|------|------|-------------|
| **Salmon** | Very fast, low memory, isoform-aware, recommended by DESeq2 authors | Transcript-level first | Most modern bulk RNA-seq |
| **STAR + featureCounts** | Gold-standard alignment, good for variant-aware or novel transcript work | High RAM (human genome ~30–40 GB), slower | When you need BAM files or junction info |

This pipeline supports both. Salmon is the default recommendation.

---

## Outputs

```
results/
├── qc/
│   └── multiqc_report.html
├── trimmed/
├── quant/                 # Salmon quant or STAR BAMs + counts
├── deseq2/
│   ├── results.csv         # full DE table
│   ├── significant.csv     # filtered by FDR + LFC
│   ├── pca.pdf
│   ├── volcano.pdf
│   ├── ma_plot.pdf
│   └── normalized_counts.csv
└── logs/
```

---

## Reference Index

You need a transcriptome (Salmon) or genome index (STAR).

**Salmon (recommended):**
```bash
# Download transcriptome (example GENCODE human)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz

salmon index -t gencode.v44.transcripts.fa.gz -i salmon_index --gencode -p 8
```

**STAR:**
```bash
STAR --runMode genomeGenerate \
  --genomeDir star_index \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile genes.gtf \
  --runThreadN 8
```

Update the path in `config/config.yaml`.

---

## Directory Structure

```
rnaseq-pipeline/
├── README.md
├── environment.yml
├── config/
│   └── config.yaml
├── data/
│   ├── samples.csv.example
│   └── raw/                 # your FASTQs go here
├── workflow/
│   └── Snakefile
├── scripts/
│   ├── run_qc.py
│   ├── run_trim.py
│   ├── run_salmon.py
│   ├── run_star.py
│   ├── deseq2_analysis.R
│   └── utils.py
├── docs/
│   └── design_and_contrasts.md
└── results/                 # generated
```

---

## Tips for Good Results

- Use biological replicates (≥3 per group recommended).
- Check the MultiQC report before proceeding to DE.
- For Salmon → DESeq2, use **tximport** (already handled in the R script).
- Include `batch` in the design formula if you have known batch effects: `~ batch + condition`.
- Always examine PCA and sample distance heatmaps for outliers.
- Report both FDR and log2 fold-change thresholds.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **Salmon** — Patro et al., Nat Methods 2017
- **DESeq2** — Love, Huber, Anders, Genome Biol 2014
- **STAR** — Dobin et al., Bioinformatics 2013
- **fastp** — Chen et al., Bioinformatics 2018
- **tximport** — Soneson et al., F1000Research 2015

---

## License

MIT
