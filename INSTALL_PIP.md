# Install notes (RNA-seq pipeline)

## Python + pip (from GitHub)

```bash
git clone https://github.com/gigi777mo/rnaseq-pipeline.git
cd rnaseq-pipeline

python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt
```

## Bioinformatics tools (not pure pip)

This pipeline also needs **Salmon or STAR**, **fastp**, **FastQC**, and R packages (**DESeq2**, **tximport**, optional **RUVSeq** / **sva**).

**Recommended:** use the included conda env for everything:

```bash
conda env create -f environment.yml
conda activate rnaseq
```

Or install tools via bioconda separately, then use pip only for the Python helpers.
