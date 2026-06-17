#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import pysam


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Make gene/transcript x cell matrices from all split BAMs in a directory."
    )
    p.add_argument("input_bam_dir", help="Directory containing split BAMs.")
    p.add_argument("original_gene_matrix_dir", help="Gene-level matrix dir with barcodes.tsv.gz and features.tsv.gz.")
    p.add_argument("original_transcript_matrix_dir", help="Transcript-level matrix dir with barcodes.tsv.gz and features.tsv.gz.")
    p.add_argument("gtf", help="GTF file used to identify chrX genes/transcripts.")

    p.add_argument("--cb-tag", default="CB")
    p.add_argument("--gn-tag", default="GN")
    p.add_argument("--tr-tag", default="TR")
    p.add_argument("--transcript-id-type", choices=["auto", "transcript_id", "transcript_name"], default="auto")

    p.add_argument("--exclude-read-ids", default=None)
    p.add_argument("--exclude-duplicates", action="store_true")
    p.add_argument("--keep-secondary", action="store_true")
    p.add_argument("--keep-supplementary", action="store_true")
    p.add_argument("--skip-transcript", action="store_true")

    p.add_argument("--outdir", default=None, help="Default: input_bam_dir/matrices")
    p.add_argument("--bam-pattern", default="*.bam")
    return p.parse_args(argv)


def load_excluded_read_ids(path: Optional[str]) -> Set[str]:
    if path is None:
        return set()

    try:
        df = pd.read_csv(path, sep="\t")
        col = "read_id" if "read_id" in df.columns else df.columns[0]
        return set(df[col].dropna().astype(str))
    except Exception:
        ids = set()
        with open(path) as f:
            for line in f:
                s = line.strip()
                if not s or s == "read_id":
                    continue
                ids.add(s.split("\t")[0])
        return ids


def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    values = {}
    for entry in str(attr).split(";"):
        entry = entry.strip()
        if not entry or " " not in entry:
            continue
        key, value = entry.split(" ", 1)
        values[key] = value.strip().strip('"')
    return values


def load_chrX_features_from_gtf(gtf_path: str) -> Tuple[Set[str], List[str], List[str]]:
    gtf_df = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"
        ],
    )

    gtf_df = gtf_df[gtf_df["seqname"].astype(str).isin(["chrX", "X"])].copy()

    gene_df = gtf_df[gtf_df["feature"] == "gene"].copy()
    gene_attrs = gene_df["attribute"].map(parse_gtf_attributes)

    chrX_genes = set(
        a.get("gene_name", a.get("gene_id"))
        for a in gene_attrs
        if a.get("gene_name", a.get("gene_id")) is not None
    )

    tx_df = gtf_df[gtf_df["feature"].isin(["transcript", "exon"])].copy()
    tx_attrs = tx_df["attribute"].map(parse_gtf_attributes)

    transcript_ids = []
    transcript_names = []
    seen_ids = set()
    seen_names = set()

    for a in tx_attrs:
        tx_id = a.get("transcript_id")
        tx_name = a.get("transcript_name")

        if tx_id and tx_id not in seen_ids:
            transcript_ids.append(tx_id)
            seen_ids.add(tx_id)

        if tx_name and tx_name not in seen_names:
            transcript_names.append(tx_name)
            seen_names.add(tx_name)

    return chrX_genes, transcript_ids, transcript_names


def read_passes_filters(read, args, excluded_read_ids: Set[str]) -> bool:
    if read.query_name in excluded_read_ids:
        return False
    if read.is_unmapped:
        return False
    if not args.keep_secondary and read.is_secondary:
        return False
    if not args.keep_supplementary and read.is_supplementary:
        return False
    if args.exclude_duplicates and read.is_duplicate:
        return False
    return True


def load_barcodes(original_matrix_dir: Path) -> List[str]:
    return pd.read_csv(original_matrix_dir / "barcodes.tsv.gz", header=None)[0].astype(str).tolist()


def load_feature_gene_order(original_matrix_dir: Path, chrX_genes: Set[str]) -> List[str]:
    features = pd.read_csv(
        original_matrix_dir / "features.tsv.gz",
        sep="\t",
        header=None
    )[1].astype(str).tolist()

    return [g for g in features if g in chrX_genes]


def format_count_matrix(
    feature_cell_counts: Dict[str, Dict[str, int]],
    row_order: List[str],
    barcodes: List[str],
) -> pd.DataFrame:

    if feature_cell_counts:
        df = pd.DataFrame.from_dict(feature_cell_counts, orient="index").fillna(0).astype(int)
    else:
        df = pd.DataFrame()

    if not df.empty:
        df.columns = [
            f"{col}-1" if not str(col).endswith("-1") else str(col)
            for col in df.columns
        ]

    observed = set(df.index.astype(str)) if not df.empty else set()
    rows = list(row_order) + sorted(observed - set(row_order))

    formatted_df = df.reindex(index=rows, columns=barcodes, fill_value=0)
    return formatted_df.fillna(0).astype(int)


def generate_matrix(
    input_bam: Path,
    barcodes: List[str],
    row_order: List[str],
    allowed_features: Optional[Set[str]],
    feature_tag: str,
    excluded_read_ids: Set[str],
    args,
) -> pd.DataFrame:

    feature_cell_counts = defaultdict(lambda: defaultdict(int))

    with pysam.AlignmentFile(str(input_bam), "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if not read_passes_filters(read, args, excluded_read_ids):
                continue

            if not read.has_tag(args.cb_tag):
                continue
            if not read.has_tag(feature_tag):
                continue

            cell_barcode = str(read.get_tag(args.cb_tag))
            feature_name = str(read.get_tag(feature_tag))

            if not feature_name or feature_name == "nan":
                continue

            if allowed_features is not None and feature_name not in allowed_features:
                continue

            feature_cell_counts[feature_name][cell_barcode] += 1

    return format_count_matrix(feature_cell_counts, row_order, barcodes)


def choose_transcript_rows(
    transcript_id_type: str,
    transcript_ids: List[str],
    transcript_names: List[str],
    bam_paths: List[Path],
    transcript_tag: str,
    args,
    excluded_read_ids: Set[str],
    max_reads: int = 250000,
) -> Tuple[List[str], Optional[Set[str]], str]:

    id_set = set(transcript_ids)
    name_set = set(transcript_names)

    if transcript_id_type == "transcript_id":
        return transcript_ids, id_set, "transcript_id"

    if transcript_id_type == "transcript_name":
        return transcript_names, name_set, "transcript_name"

    observed = set()
    n_seen = 0

    for bam_path in bam_paths:
        with pysam.AlignmentFile(str(bam_path), "rb") as bamfile:
            for read in bamfile.fetch(until_eof=True):
                if not read_passes_filters(read, args, excluded_read_ids):
                    continue

                if read.has_tag(transcript_tag):
                    observed.add(str(read.get_tag(transcript_tag)))

                n_seen += 1
                if n_seen >= max_reads:
                    break

        if n_seen >= max_reads:
            break

    id_hits = len(observed & id_set)
    name_hits = len(observed & name_set)

    if id_hits >= name_hits and id_hits > 0:
        return transcript_ids, id_set, "transcript_id"

    if name_hits > 0:
        return transcript_names, name_set, "transcript_name"

    return [], None, "observed_transcript_tag_values"


def matrix_name_from_bam(bam_path: Path, transcript: bool = False) -> str:
    name = bam_path.name

    if name.endswith(".bam"):
        name = name[:-4]

    suffix = name.split(".")[-1]
    known = {"X1", "X2", "Xa", "Xi", "low", "amb", "unknown"}

    if suffix in known:
        base = suffix
    else:
        base = name

    if transcript:
        return f"transcript_{base}.csv"

    return f"{base}.csv"


def main(argv=None):
    args = parse_args(argv)

    input_bam_dir = Path(args.input_bam_dir)
    original_gene_matrix_dir = Path(args.original_gene_matrix_dir)
    original_transcript_matrix_dir = Path(args.original_transcript_matrix_dir)

    outdir = Path(args.outdir) if args.outdir else input_bam_dir / "matrices"
    outdir.mkdir(parents=True, exist_ok=True)

    bam_paths = sorted(input_bam_dir.glob(args.bam_pattern))
    bam_paths = [p for p in bam_paths if p.is_file() and p.suffix == ".bam"]

    if not bam_paths:
        raise FileNotFoundError(f"No BAM files found in {input_bam_dir} using pattern {args.bam_pattern}")

    excluded_read_ids = load_excluded_read_ids(args.exclude_read_ids)
    print(f"Loaded {len(excluded_read_ids):,} excluded read IDs")

    chrX_genes, chrX_transcript_ids, chrX_transcript_names = load_chrX_features_from_gtf(args.gtf)
    print(f"Loaded {len(chrX_genes):,} chrX genes from GTF")
    print(f"Loaded {len(chrX_transcript_ids):,} chrX transcript IDs from GTF")
    print(f"Loaded {len(chrX_transcript_names):,} chrX transcript names from GTF")

    print(f"Found {len(bam_paths)} BAMs:")
    for bam_path in bam_paths:
        print(f"  {bam_path}")

    gene_barcodes = load_barcodes(original_gene_matrix_dir)
    gene_rows = load_feature_gene_order(original_gene_matrix_dir, chrX_genes)

    print("\nMaking gene-level matrices")
    for bam_path in bam_paths:
        out_csv = outdir / matrix_name_from_bam(bam_path, transcript=False)

        print(f"Counting gene-level BAM: {bam_path}")
        matrix = generate_matrix(
            bam_path,
            gene_barcodes,
            gene_rows,
            chrX_genes,
            args.gn_tag,
            excluded_read_ids,
            args,
        )

        matrix.to_csv(out_csv)
        print(f"  wrote {out_csv} shape={matrix.shape}")

    if not args.skip_transcript:
        transcript_barcodes = load_barcodes(original_transcript_matrix_dir)

        tx_rows, tx_allowed, tx_mode = choose_transcript_rows(
            args.transcript_id_type,
            chrX_transcript_ids,
            chrX_transcript_names,
            bam_paths,
            args.tr_tag,
            args,
            excluded_read_ids,
        )

        print(f"\nTranscript row mode: {tx_mode}")

        if tx_allowed is None:
            print(
                "Warning: transcript BAM tag values did not match GTF transcript_id "
                "or transcript_name in sampled reads. Transcript matrices will use "
                "observed transcript tag values only, without chrX restriction."
            )

        print("\nMaking transcript-level matrices")
        for bam_path in bam_paths:
            out_csv = outdir / matrix_name_from_bam(bam_path, transcript=True)

            print(f"Counting transcript-level BAM: {bam_path}")
            matrix = generate_matrix(
                bam_path,
                transcript_barcodes,
                tx_rows,
                tx_allowed,
                args.tr_tag,
                excluded_read_ids,
                args,
            )

            matrix.to_csv(out_csv)
            print(f"  wrote {out_csv} shape={matrix.shape}")

    print(f"\nDone. Output directory: {outdir}")


if __name__ == "__main__":
    raise SystemExit(main())