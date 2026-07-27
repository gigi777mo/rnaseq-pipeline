# START HERE (no bioinformatics experience needed)

---

> ## 🔴 USE MINICONDA — REQUIRED
>
> **You must use [Miniconda](https://docs.conda.io/en/latest/miniconda.html).**  
> This pipeline needs Salmon/STAR, fastp, and R packages.  
> **Pip alone will not install those tools.**  
> Skip Miniconda → expect errors.

---

## What this does (plain English)

Takes RNA-seq **FASTQ** files and:

1. Checks quality  
2. Trims adapters  
3. Counts genes (Salmon)  
4. Finds which genes change between conditions (DESeq2)  
5. Optional: spike-ins, RUVSeq, SVA batch correction  

---

## Step 1 — Install Miniconda (one time)

https://docs.conda.io/en/latest/miniconda.html  
Install → **close and reopen** the terminal / Anaconda Prompt.

---

## Step 2 — Get the code

```bash
git clone https://github.com/gigi777mo/rnaseq-pipeline.git
cd rnaseq-pipeline
```

Or Download ZIP from GitHub.

---

## Step 3 — Create environment (one time, can take a while)

```bash
conda env create -f environment.yml
conda activate rnaseq
```

Type `y` if asked. Wait until it finishes.

---

## Step 4 — Your data

1. Put FASTQ files in `data/raw/`  
2. Edit `data/samples.csv` (copy from `data/samples.csv.example`)  
3. Edit `config/config.yaml` (genome/transcript index paths if needed)  

Ask someone in the lab for the **transcriptome index** path if you do not have one yet.

---

## Step 5 — Run

```bash
conda activate rnaseq
snakemake --cores 4 --use-conda
```

Or run individual Python/R scripts listed in the README when you only need one step.

---

## Step 6 — Results

Look under `results/` for QC reports, quantifications, and DESeq2 tables (paths depend on your config).

---

## If it fails

| Problem | Fix |
|---------|-----|
| `conda not found` | Install Miniconda; **new** terminal |
| Env create fails | Update conda: `conda update -n base conda` and retry |
| Missing index | Set correct Salmon/STAR index in `config/config.yaml` |
| Out of memory | Lower `--cores` or run on a bigger machine |

**Do not switch to pip-only for this pipeline.**  
Citations: [docs/CITATIONS.md](docs/CITATIONS.md)
