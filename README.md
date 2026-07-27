# RNA-seq Pipeline

---

> ## 🔴 USE MINICONDA — REQUIRED
>
> **Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.**  
> `conda env create -f environment.yml` → `conda activate rnaseq`  
> **Pip cannot replace Salmon/STAR/DESeq2 here.**

---

**New user?** Open **[START_HERE.md](START_HERE.md)** and follow the steps only.

Bulk RNA-seq: QC → trim → Salmon/STAR → DESeq2, with optional spike-ins, RUVSeq, and SVA.

```bash
conda env create -f environment.yml
conda activate rnaseq
snakemake --cores 4 --use-conda
```

Citations: [docs/CITATIONS.md](docs/CITATIONS.md)

## License

MIT
