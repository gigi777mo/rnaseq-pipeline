#!/usr/bin/env python3
"""Trim paired-end FASTQs with fastp."""

import argparse
import subprocess
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", default="results/trimmed")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    indir = Path(args.input)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # Pair R1/R2 by common prefix
    r1_files = sorted(indir.glob("*_R1*.fastq.gz")) + sorted(indir.glob("*_1*.fastq.gz"))
    for r1 in r1_files:
        name = r1.name.replace("_R1.fastq.gz", "").replace("_1.fastq.gz", "").replace("_R1_001.fastq.gz", "")
        # Try common pairing patterns
        candidates = [
            indir / f"{name}_R2.fastq.gz",
            indir / f"{name}_2.fastq.gz",
            indir / f"{name}_R2_001.fastq.gz",
        ]
        r2 = next((c for c in candidates if c.exists()), None)
        if r2 is None:
            print(f"[!] No R2 found for {r1.name}, skipping")
            continue

        out_r1 = outdir / f"{name}_R1.fastq.gz"
        out_r2 = outdir / f"{name}_R2.fastq.gz"
        html = outdir / f"{name}.fastp.html"
        json = outdir / f"{name}.fastp.json"

        print(f"[+] Trimming {name}...")
        cmd = [
            "fastp",
            "-i", str(r1), "-I", str(r2),
            "-o", str(out_r1), "-O", str(out_r2),
            "-h", str(html), "-j", str(json),
            "--detect_adapter_for_pe",
            "-w", str(args.threads)
        ]
        subprocess.run(cmd, check=True)

    print("[+] Trimming complete")

if __name__ == "__main__":
    main()
