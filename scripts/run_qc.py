#!/usr/bin/env python3
"""
Enhanced QC for RNA-seq:
- FastQC on raw and/or trimmed reads
- MultiQC aggregation
- Basic read-count warnings
- Optional spike-in metrics (delegates to R script when quantifications exist)
"""

import argparse
import subprocess
import sys
from pathlib import Path

def run_fastqc(fastqs, outdir, threads):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not fastqs:
        print("[!] No FASTQ files provided to FastQC")
        return
    print(f"[+] FastQC on {len(fastqs)} files → {outdir}")
    subprocess.run(
        ["fastqc", "-t", str(threads), "-o", str(outdir)] + [str(f) for f in fastqs],
        check=True
    )

def run_multiqc(qc_dirs, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["multiqc"] + [str(d) for d in qc_dirs] + ["-o", str(outdir), "--force"]
    print("[+] MultiQC...")
    subprocess.run(cmd, check=True)
    print(f"[+] MultiQC report: {outdir}/multiqc_report.html")

def count_reads_fastq(path):
    """Approximate read count from a gzipped FASTQ (counts lines / 4)."""
    import gzip
    n = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for _ in f:
            n += 1
    return n // 4

def main():
    parser = argparse.ArgumentParser(description="RNA-seq QC (FastQC + MultiQC + basic checks)")
    parser.add_argument("--input", required=True, help="Directory with FASTQ files (raw or trimmed)")
    parser.add_argument("--out", default="results/qc")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-reads", type=int, default=1_000_000,
                        help="Warn if a sample has fewer reads than this")
    parser.add_argument("--stage", default="raw", choices=["raw", "trimmed", "both"],
                        help="Which stage the input FASTQs belong to")
    parser.add_argument("--salmon-dir", default=None,
                        help="If provided and spike-ins were used, run spike-in QC")
    parser.add_argument("--samples", default="data/samples.csv")
    parser.add_argument("--spikein-prefix", default="ERCC-")
    args = parser.parse_args()

    indir = Path(args.input)
    outdir = Path(args.out)
    stage_dir = outdir / args.stage
    stage_dir.mkdir(parents=True, exist_ok=True)

    fastqs = sorted(list(indir.glob("*.fastq.gz")) + list(indir.glob("*.fq.gz")))
    if not fastqs:
        sys.exit(f"No FASTQ files found in {indir}")

    # FastQC
    run_fastqc(fastqs, stage_dir, args.threads)

    # Basic read-count check (on R1 only to avoid double-counting)
    r1s = [f for f in fastqs if "_R1" in f.name or "_1" in f.name]
    print("\n[+] Approximate read counts (R1):")
    for f in r1s:
        try:
            n = count_reads_fastq(f)
            status = "OK" if n >= args.min_reads else "LOW"
            print(f"  {f.name}: {n:,} reads  [{status}]")
            if n < args.min_reads:
                print(f"    ⚠️  Below --min-reads threshold ({args.min_reads:,})")
        except Exception as e:
            print(f"  {f.name}: could not count ({e})")

    # MultiQC
    run_multiqc([stage_dir], outdir)

    # Optional spike-in QC
    if args.salmon_dir and Path(args.salmon_dir).exists():
        print("[+] Running spike-in QC...")
        cmd = [
            "Rscript", "scripts/spikein_qc.R",
            "--salmon_dir", args.salmon_dir,
            "--samples", args.samples,
            "--prefix", args.spikein_prefix,
            "--out", str(outdir)
        ]
        subprocess.run(cmd, check=False)

    print("\n[+] QC finished.")

if __name__ == "__main__":
    main()
