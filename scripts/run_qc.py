#!/usr/bin/env python3
"""Run FastQC + MultiQC on raw FASTQ files."""

import argparse
import subprocess
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Directory with raw FASTQ files")
    parser.add_argument("--out", default="results/qc")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    indir = Path(args.input)
    outdir = Path(args.out)
    raw_qc = outdir / "raw"
    raw_qc.mkdir(parents=True, exist_ok=True)

    fastqs = list(indir.glob("*.fastq.gz")) + list(indir.glob("*.fq.gz"))
    if not fastqs:
        raise SystemExit(f"No FASTQ files found in {indir}")

    print(f"[+] Running FastQC on {len(fastqs)} files...")
    subprocess.run(["fastqc", "-t", str(args.threads), "-o", str(raw_qc)] + [str(f) for f in fastqs], check=True)

    print("[+] Running MultiQC...")
    subprocess.run(["multiqc", str(raw_qc), "-o", str(outdir), "--force"], check=True)
    print(f"[+] MultiQC report: {outdir}/multiqc_report.html")

if __name__ == "__main__":
    main()
