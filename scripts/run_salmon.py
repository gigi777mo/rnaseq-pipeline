#!/usr/bin/env python3
"""Quantify transcripts with Salmon."""

import argparse
import subprocess
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Directory with trimmed FASTQs")
    parser.add_argument("--index", required=True, help="Salmon index directory")
    parser.add_argument("--out", default="results/salmon")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--libtype", default="A")
    args = parser.parse_args()

    indir = Path(args.input)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    r1_files = sorted(indir.glob("*_R1.fastq.gz"))
    for r1 in r1_files:
        name = r1.name.replace("_R1.fastq.gz", "")
        r2 = indir / f"{name}_R2.fastq.gz"
        if not r2.exists():
            print(f"[!] Missing R2 for {name}")
            continue

        sample_out = outdir / name
        print(f"[+] Salmon quant: {name}")
        cmd = [
            "salmon", "quant",
            "-i", args.index,
            "-l", args.libtype,
            "-1", str(r1), "-2", str(r2),
            "-o", str(sample_out),
            "-p", str(args.threads),
            "--gcBias", "--seqBias"
        ]
        subprocess.run(cmd, check=True)

    print("[+] Salmon quantification complete")

if __name__ == "__main__":
    main()
