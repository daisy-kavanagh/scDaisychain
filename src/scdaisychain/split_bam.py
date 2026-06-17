#!/usr/bin/env python3
"""Split a tagged BAM into X1/X2/Xa/Xi/low/amb/unknown BAMs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pysam


LABELS = ["X1", "X2", "Xa", "Xi", "low", "amb", "unknown"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Split a BAM tagged by xci-tag-bam into X1/X2/Xa/Xi/low/amb/unknown BAMs."
    )
    parser.add_argument("inbam", help="Input tagged BAM/CRAM.")
    parser.add_argument("outdir", help="Output directory.")
    parser.add_argument("--prefix", default="split", help="Output BAM prefix. Default: split.")
    return parser.parse_args(argv)


def split_bam(args) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bam = pysam.AlignmentFile(args.inbam, "rb")

    paths = {k: outdir / f"{args.prefix}.{k}.bam" for k in LABELS}
    outfiles = {k: pysam.AlignmentFile(str(path), "wb", header=bam.header) for k, path in paths.items()}
    counts = {k: 0 for k in LABELS}

    try:
        for read in bam.fetch(until_eof=True):
            try:
                ma = read.get_tag("MA")  # H1/H2/low/amb/none
            except KeyError:
                ma = "unknown"

            try:
                xx = read.get_tag("XX")  # Xa/Xi/low/amb/no_CB/no_active_X/etc
            except KeyError:
                xx = "unknown"

            # Split by haplotype X1 / X2.
            if ma == "H1":
                outfiles["X1"].write(read)
                counts["X1"] += 1
            elif ma == "H2":
                outfiles["X2"].write(read)
                counts["X2"] += 1

            # Split by Xa / Xi / low / amb / unknown.
            if xx == "Xa":
                outfiles["Xa"].write(read)
                counts["Xa"] += 1
            elif xx == "Xi":
                outfiles["Xi"].write(read)
                counts["Xi"] += 1
            elif xx == "low":
                outfiles["low"].write(read)
                counts["low"] += 1
            elif xx == "amb":
                outfiles["amb"].write(read)
                counts["amb"] += 1
            else:
                outfiles["unknown"].write(read)
                counts["unknown"] += 1
    finally:
        bam.close()
        for f in outfiles.values():
            f.close()

    for name, path in paths.items():
        pysam.index(str(path))
        print(f"{name}\t{counts[name]}\t{path}\t{path}.bai")


def main(argv=None):
    args = parse_args(argv)
    split_bam(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
