#!/usr/bin/env python3
"""Tag BAM reads with haplotype (MA/HP) and Xa/Xi status (XX)."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Tuple

import pysam


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Tag reads using a phased VCF and per-cell active-X table. "
            "Adds MS/NS/MC/NC/MI, MA, optional HP, and XX tags."
        )
    )
    parser.add_argument("bam", help="Input BAM/CRAM.")
    parser.add_argument("vcf", help="Phased VCF/BCF. Must contain a phased GT sample.")
    parser.add_argument("active_x_csv", help="CSV with columns cell_barcode and active_X.")
    parser.add_argument("outbam", help="Output tagged BAM.")
    parser.add_argument(
        "--mode",
        choices=["quality", "count"],
        default="quality",
        help="Assign haplotype by summed base quality or SNP count. Default: quality.",
    )
    parser.add_argument("--min-base-qual", type=int, default=5)
    parser.add_argument("--min-score-diff", type=int, default=10)
    parser.add_argument("--min-total-score", type=int, default=10)
    parser.add_argument("--min-count-diff", type=int, default=1)
    parser.add_argument("--min-total-count", type=int, default=1)
    return parser.parse_args(argv)


def load_active_x(path: str) -> Dict[str, str]:
    active_x: Dict[str, str] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"cell_barcode", "active_X"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        for row in reader:
            active_x[row["cell_barcode"]] = row["active_X"]
    return active_x


def load_phased_snps(vcf_path: str) -> Dict[Tuple[str, int], Tuple[str, str]]:
    vcf = pysam.VariantFile(vcf_path)
    try:
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError(f"VCF has no samples: {vcf_path}")
        sample = samples[0]

        snp_dict: Dict[Tuple[str, int], Tuple[str, str]] = {}
        for rec in vcf.fetch():
            if len(rec.ref) != 1 or len(rec.alts or []) != 1 or len(rec.alts[0]) != 1:
                continue

            gt = rec.samples[sample].get("GT")
            if gt is None or len(gt) != 2:
                continue
            if not rec.samples[sample].phased:
                continue
            if gt[0] is None or gt[1] is None:
                continue

            alleles = [rec.ref.upper(), rec.alts[0].upper()]
            try:
                snp_dict[(rec.chrom, rec.pos)] = (alleles[gt[0]], alleles[gt[1]])
            except IndexError:
                continue
        return snp_dict
    finally:
        vcf.close()


def tag_bam(args) -> None:
    active_x = load_active_x(args.active_x_csv)
    snp_dict = load_phased_snps(args.vcf)

    out_path = Path(args.outbam)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    bam = pysam.AlignmentFile(args.bam, "rb")
    out = pysam.AlignmentFile(args.outbam, "wb", header=bam.header)

    try:
        for read in bam.fetch(until_eof=True):
            h1_score = h2_score = 0
            h1_count = h2_count = informative = 0

            if not read.is_unmapped:
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if query_pos is None or ref_pos is None:
                        continue
                    key = (read.reference_name, ref_pos + 1)

                    if key not in snp_dict:
                        continue
                    if read.query_sequence is None or read.query_qualities is None:
                        continue

                    base = read.query_sequence[query_pos].upper()
                    qual = read.query_qualities[query_pos]

                    if qual < args.min_base_qual:
                        continue

                    h1, h2 = snp_dict[key]

                    if base == h1:
                        h1_score += qual
                        h1_count += 1
                        informative += 1
                    elif base == h2:
                        h2_score += qual
                        h2_count += 1
                        informative += 1

            read.set_tag("MS", int(h1_score))     # H1 quality score
            read.set_tag("NS", int(h2_score))     # H2 quality score
            read.set_tag("MC", int(h1_count))     # H1 SNP count
            read.set_tag("NC", int(h2_count))     # H2 SNP count
            read.set_tag("MI", int(informative))  # informative SNP count

            hap_assignment = "none"

            if informative == 0:
                read.set_tag("MA", "none")

            elif args.mode == "quality":
                total_score = h1_score + h2_score
                score_diff = abs(h1_score - h2_score)

                if total_score < args.min_total_score:
                    read.set_tag("MA", "low")
                    hap_assignment = "low"
                elif score_diff < args.min_score_diff:
                    read.set_tag("MA", "amb")
                    hap_assignment = "amb"
                elif h1_score > h2_score:
                    read.set_tag("MA", "H1")
                    read.set_tag("HP", 1)
                    hap_assignment = "H1"
                else:
                    read.set_tag("MA", "H2")
                    read.set_tag("HP", 2)
                    hap_assignment = "H2"

            elif args.mode == "count":
                total_count = h1_count + h2_count
                count_diff = abs(h1_count - h2_count)

                if total_count < args.min_total_count:
                    read.set_tag("MA", "low")
                    hap_assignment = "low"
                elif count_diff < args.min_count_diff:
                    read.set_tag("MA", "amb")
                    hap_assignment = "amb"
                elif h1_count > h2_count:
                    read.set_tag("MA", "H1")
                    read.set_tag("HP", 1)
                    hap_assignment = "H1"
                else:
                    read.set_tag("MA", "H2")
                    read.set_tag("HP", 2)
                    hap_assignment = "H2"

            x_status = "none"

            try:
                cb = read.get_tag("CB")
            except KeyError:
                cb = None

            if hap_assignment in {"H1", "H2"} and cb in active_x:
                ax = active_x[cb]

                if ax == "X1":
                    x_status = "Xa" if hap_assignment == "H1" else "Xi"
                elif ax == "X2":
                    x_status = "Xa" if hap_assignment == "H2" else "Xi"
                else:
                    x_status = "unknown"

            elif hap_assignment in {"low", "amb"}:
                x_status = hap_assignment
            elif cb is None:
                x_status = "no_CB"
            elif cb not in active_x:
                x_status = "no_active_X"

            read.set_tag("XX", x_status)
            out.write(read)
    finally:
        bam.close()
        out.close()


def main(argv=None):
    args = parse_args(argv)
    tag_bam(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
