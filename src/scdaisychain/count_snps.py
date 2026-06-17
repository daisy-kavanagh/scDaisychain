#!/usr/bin/env python3
"""
Per-cell ref/alt/other counts at SNPs from a VCF, using CB (cell barcode)
and GN (gene) tags in BAM/CRAM reads.

Output count columns:
  cell_barcode, contig, position, gene, refAllele, altAllele,
  refCount, altCount, totalCount, otherCount

Optional GTF annotation/filtering adds columns:
  gtf_gene, gtf_gene_ids, gtf_gene_n, gtf_multi_gene, gtf_match_type,
  gtf_min_distance, gtf_gene_conflict, gtf_conflict_reason,
  tsv_gene_in_gtf_list

Optional dropped-read outputs:
  --dropped-read-ids-out          read x SNP evidence rows for dropped rows
  --unique-dropped-read-ids-out   one unique read ID per line
  --dropped-read-summary-out      one row per dropped read with dropped SNP/gene summaries

Main implementation notes:
- One BAM pileup pass per region/window, rather than one pileup call per SNP.
- Deterministic parallel merge using future metadata.
- Temp files are cleaned up even if a worker/merge step fails.
"""

import argparse
import os
import sys
import tempfile
from collections import defaultdict
import gzip
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple
import bisect

import pysam
import pandas as pd


HEADER_COLS = [
    "cell_barcode", "contig", "position", "gene", "refAllele", "altAllele",
    "refCount", "altCount", "totalCount", "otherCount",
]

EVIDENCE_COLS = [
    "read_id", "cell_barcode", "contig", "position", "gene",
    "refAllele", "altAllele", "allele_class",
]


# --------------------------- CLI ---------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=(
            "Per-cell ref/alt/other counts at SNPs from VCF using CB (cell barcode) "
            "and GN (gene) tags in BAM/CRAM reads, with optional GTF filtering."
        )
    )

    p.add_argument("--bam", required=True, help="Coordinate-sorted, indexed BAM/CRAM (.bai/.crai required).")
    p.add_argument("--vcf", required=True, help="VCF/BCF of SNPs (indexed .tbi/.csi required for region fetch).")
    p.add_argument("--out", required=True, help="Output TSV path. If --gtf is supplied, this is the annotated/final filtered count TSV.")

    p.add_argument("--gtf", default=None, help="Optional GTF file for gene-position annotation and filtering.")
    p.add_argument("--slop", type=int, default=1000, help="Generous bp padding around GTF gene intervals. Default: 1000.")
    p.add_argument("--max-nearest-dist", type=int, default=10000, help="If no overlap/slop hit, assign nearest gene up to this distance. Default: 10000. Use 0 to disable.")
    p.add_argument("--gtf-feature", default="gene", help="GTF feature to use. Default: gene.")
    p.add_argument("--gene-name-attr", default="gene_name", help="GTF attribute for gene symbol. Default: gene_name.")
    p.add_argument("--gene-id-attr", default="gene_id", help="GTF attribute for gene ID. Default: gene_id.")

    p.add_argument("--drop-conflicts", action="store_true", help="After annotation, drop rows where gtf_gene_conflict is True.")
    p.add_argument(
        "--drop-multi-tsv-and-gtf",
        action="store_true",
        help="After optional conflict removal, drop entire SNP positions that still have 2+ nonblank TSV genes and 2+ GTF genes.",
    )
    p.add_argument("--summary-out", default=None, help="Optional summary TSV path. Default: <out>.summary.tsv")

    p.add_argument(
        "--dropped-read-ids-out",
        default=None,
        help="Optional TSV of read x SNP evidence rows supporting rows removed by filters. Requires --gtf and at least one drop option.",
    )
    p.add_argument(
        "--unique-dropped-read-ids-out",
        default=None,
        help="Optional file containing one unique dropped read_id per line. Requires dropped-read evidence generation.",
    )
    p.add_argument(
        "--dropped-read-summary-out",
        default=None,
        help="Optional TSV with one row per dropped read: n_dropped_snps, cells, genes, drop_reasons, allele_classes.",
    )
    p.add_argument(
        "--read-evidence-out",
        default=None,
        help="Optional TSV of all read-level SNP evidence before aggregation. Can be large.",
    )

    p.add_argument("--cb-tag", default="CB", help="Read tag for cell barcode (default: CB).")
    p.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality to count a read (default: 20).")
    p.add_argument("--min-bq", type=int, default=10, help="Minimum base quality to count a base (default: 13).")
    p.add_argument("--max-depth", type=int, default=100000, help="Max pileup depth per site (default: 100000).")
    p.add_argument("--allow-secondary", action="store_true", help="Include secondary/supplementary alignments.")
    p.add_argument("--count-duplicates", action="store_true", help="Include reads marked as duplicates.")
    p.add_argument("--skip-indels", action="store_true", help="Skip positions where read has an indel/deletion at the site.")
    p.add_argument("--ignore-overlaps", action="store_true", help="Pass ignore_overlaps=True to pileup.")

    p.add_argument("--region", default=None, help="Region like 'chrX' or 'chrX:1-20000000'.")
    p.add_argument("--chrom", default=None, help="Chromosome to process/shard, e.g. chrX.")
    p.add_argument("--window", type=int, default=5_000_000, help="Window size for region sharding (default: 5 Mb).")

    p.add_argument("--procs", type=int, default=1, help="Parallel region workers (default: 1).")
    p.add_argument("--threads", type=int, default=1, help="HTSlib BGZF threads per worker (default: 1).")

    p.add_argument("--het-only", action="store_true", help="Restrict to heterozygous SNPs (requires genotype data).")
    p.add_argument("--sample", default=None, help="Sample name for GT. Default: first sample in the VCF.")

    p.add_argument("--debug", action="store_true", help="Enable debug prints.")
    p.add_argument("--debug-examples", type=int, default=5, help="Max example reads to print per region (default: 5).")

    args = p.parse_args(argv)

    if args.window <= 0:
        p.error("--window must be > 0")
    if args.procs <= 0:
        p.error("--procs must be > 0")
    if args.threads <= 0:
        p.error("--threads must be > 0")
    if args.region and args.chrom:
        p.error("Use only one of --region or --chrom.")
    if args.procs > 1 and not (args.region or args.chrom):
        p.error("With --procs > 1, provide --region or --chrom so work can be safely sharded.")

    wants_dropped_outputs = bool(args.dropped_read_ids_out or args.unique_dropped_read_ids_out or args.dropped_read_summary_out)
    if (args.drop_conflicts or args.drop_multi_tsv_and_gtf or wants_dropped_outputs) and not args.gtf:
        p.error("--gtf is required when using filtering or dropped-read outputs")
    if wants_dropped_outputs and not (args.drop_conflicts or args.drop_multi_tsv_and_gtf):
        p.error("Dropped-read outputs only make sense with --drop-conflicts and/or --drop-multi-tsv-and-gtf")

    return args


# ------------------------ Utilities ------------------------

def is_snp(rec) -> bool:
    if rec.alts is None or len(rec.alts) == 0:
        return False
    return len(rec.ref) == 1 and len(rec.alts[0]) == 1


def parse_region(region_str: str) -> Tuple[str, int, Optional[int]]:
    if ":" not in region_str:
        return region_str, 1, None
    chrom, rest = region_str.split(":", 1)
    rest = rest.replace(",", "")
    if "-" not in rest:
        return chrom, int(rest), None
    s, e = rest.split("-", 1)
    return chrom, int(s), int(e)


def region_to_fetch_args(region: str) -> Tuple[str, int, int]:
    chrom, start1, end1 = parse_region(region)
    if end1 is None:
        raise ValueError(f"Internal error: region {region!r} does not have an end coordinate")
    return chrom, start1 - 1, end1


def chrom_length_from_bam(bam_path: str, chrom: str) -> int:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        try:
            return bam.get_reference_length(chrom)
        except (KeyError, ValueError):
            pass
        for sq in bam.header.to_dict().get("SQ", []):
            if sq.get("SN") == chrom:
                return int(sq["LN"])
    raise ValueError(f"Chromosome {chrom!r} not found in BAM header.")


def bam_contig_lengths(bam_path: str) -> Dict[str, int]:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return dict(zip(bam.references, bam.lengths))


def vcf_contigs(vcf_path: str) -> List[str]:
    with pysam.VariantFile(vcf_path) as vcf:
        return list(vcf.header.contigs)


def make_windows(chrom: str, length: int, win: int, start: int = 1, end: Optional[int] = None) -> List[str]:
    last = length if end is None else min(end, length)
    out = []
    s = max(1, start)
    while s <= last:
        e = min(s + win - 1, last)
        out.append(f"{chrom}:{s}-{e}")
        s = e + 1
    return out


def build_flag_filter(allow_secondary: bool, count_duplicates: bool) -> int:
    mask = 0
    if not allow_secondary:
        mask |= pysam.FSECONDARY | pysam.FSUPPLEMENTARY
    if not count_duplicates:
        mask |= pysam.FDUP
    return mask


def genotype_is_het_01(rec, sample_name: Optional[str]) -> bool:
    if sample_name is None:
        return False
    smp = rec.samples.get(sample_name)
    if smp is None or "GT" not in smp:
        return False
    gt = smp["GT"]
    if not gt or len(gt) < 2:
        return False
    try:
        return set(gt) == {0, 1}
    except TypeError:
        return False


def get_sample_for_het_only(vcf, requested: Optional[str]) -> Optional[str]:
    if requested is not None:
        if requested not in vcf.header.samples:
            raise ValueError(f"Requested --sample {requested!r} not found in VCF samples: {list(vcf.header.samples)}")
        return requested
    if len(vcf.header.samples) == 0:
        return None
    return list(vcf.header.samples)[0]


def build_regions(args) -> List[str]:
    if args.region:
        chrom, start, end = parse_region(args.region)
        L = chrom_length_from_bam(args.bam, chrom)
        if end is None:
            end = L
        if args.procs > 1:
            return make_windows(chrom, L, args.window, start=start, end=end)
        return [f"{chrom}:{start}-{min(end, L)}"]

    if args.chrom:
        L = chrom_length_from_bam(args.bam, args.chrom)
        if args.procs > 1:
            return make_windows(args.chrom, L, args.window)
        return [f"{args.chrom}:1-{L}"]

    lengths = bam_contig_lengths(args.bam)
    contigs = vcf_contigs(args.vcf)
    if not contigs:
        raise ValueError("VCF header has no contig lines. Please provide --chrom or --region.")
    regions = [f"{c}:1-{lengths[c]}" for c in contigs if c in lengths]
    if not regions:
        raise ValueError("No VCF contigs were found in the BAM header. Check chr naming or provide --chrom/--region.")
    return regions


# --------------------- GTF annotation / filtering ---------------------

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def parse_gtf_attrs(attr_string: str) -> Dict[str, str]:
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_string):
        out[m.group(1)] = m.group(2)
    return out


def norm_gene(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if s in {"", ".", "-"} or s.upper() in {"NA", "NAN", "NONE"}:
        return ""
    return s


def load_gtf_genes(gtf_path: str, feature: str, gene_name_attr: str, gene_id_attr: str):
    genes_by_contig = defaultdict(list)
    n = 0
    with open_text(gtf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig, source, feat, start, end, score, strand, frame, attrs = parts
            if feat != feature:
                continue
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            a = parse_gtf_attrs(attrs)
            gene_id = a.get(gene_id_attr, "")
            gene_name = a.get(gene_name_attr) or a.get("gene") or gene_id
            if not gene_name and not gene_id:
                continue
            genes_by_contig[contig].append({
                "contig": contig,
                "start": start,
                "end": end,
                "gene_name": gene_name,
                "gene_id": gene_id,
                "strand": strand,
            })
            n += 1

    index = {}
    for contig, rows in genes_by_contig.items():
        rows = sorted(rows, key=lambda r: (r["start"], r["end"], r["gene_name"]))
        index[contig] = {"rows": rows, "starts": [r["start"] for r in rows]}
    print(f"Loaded {n:,} {feature!r} features from GTF across {len(index):,} contigs.", file=sys.stderr)
    return index


def candidate_contigs(contig: str):
    c = str(contig)
    out = [c]
    out.append(c[3:] if c.startswith("chr") else "chr" + c)
    if c == "MT":
        out += ["chrM", "M"]
    if c == "chrM":
        out += ["MT", "M"]
    seen = set()
    uniq = []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq


def distance_to_interval(pos: int, start: int, end: int) -> int:
    if start <= pos <= end:
        return 0
    return start - pos if pos < start else pos - end


def annotate_pos(contig: str, pos: int, gtf_index, slop: int, max_nearest_dist: int):
    data = None
    for c in candidate_contigs(contig):
        if c in gtf_index:
            data = gtf_index[c]
            break
    if data is None:
        return {
            "gtf_gene": "", "gtf_gene_ids": "", "gtf_gene_n": 0,
            "gtf_multi_gene": False, "gtf_match_type": "none", "gtf_min_distance": pd.NA,
        }

    rows = data["rows"]
    starts = data["starts"]
    right = bisect.bisect_right(starts, pos + slop)
    hits = []
    exact = []

    for r in rows[:right]:
        if r["end"] < pos - slop:
            continue
        d = distance_to_interval(pos, r["start"], r["end"])
        if d == 0:
            exact.append((d, r))
        if d <= slop:
            hits.append((d, r))

    if exact:
        assigned = exact
        match_type = "exact_overlap"
    elif hits:
        assigned = hits
        match_type = "within_slop"
    else:
        assigned = []
        match_type = "none"
        if max_nearest_dist and max_nearest_dist > 0:
            idx = bisect.bisect_left(starts, pos)
            lo = max(0, idx - 200)
            hi = min(len(rows), idx + 200)
            nearest = []
            best_d = None
            for r in rows[lo:hi]:
                d = distance_to_interval(pos, r["start"], r["end"])
                if best_d is None or d < best_d:
                    best_d = d
                    nearest = [(d, r)]
                elif d == best_d:
                    nearest.append((d, r))
            if best_d is not None and best_d <= max_nearest_dist:
                assigned = nearest
                match_type = "nearest_within_max"

    seen = set()
    dedup = []
    for d, r in sorted(assigned, key=lambda x: (x[0], x[1]["start"], x[1]["end"], x[1]["gene_name"])):
        key = (r["gene_id"], r["gene_name"])
        if key in seen:
            continue
        seen.add(key)
        dedup.append((d, r))

    genes = [r["gene_name"] for d, r in dedup]
    gene_ids = [r["gene_id"] for d, r in dedup]
    min_d = min([d for d, r in dedup], default=pd.NA)
    return {
        "gtf_gene": ";".join(genes),
        "gtf_gene_ids": ";".join(gene_ids),
        "gtf_gene_n": len(genes),
        "gtf_multi_gene": len(genes) > 1,
        "gtf_match_type": match_type,
        "gtf_min_distance": min_d,
    }


def conflict_reason(tsv_gene: str, ann: Dict) -> Tuple[bool, str, bool]:
    tsv_gene = norm_gene(tsv_gene)
    gtf_genes = [g for g in str(ann.get("gtf_gene", "")).split(";") if g]
    if not tsv_gene:
        return (False, "tsv_gene_blank", False) if gtf_genes else (False, "both_unannotated", False)
    if not gtf_genes:
        return True, "no_gtf_gene_near_position", False
    in_list = tsv_gene in set(gtf_genes)
    if in_list:
        return False, "tsv_gene_matches_gtf", True
    if ann.get("gtf_match_type") in {"exact_overlap", "within_slop"}:
        return True, "tsv_gene_not_in_overlapping_or_slop_gtf_genes", False
    if ann.get("gtf_match_type") == "nearest_within_max":
        return True, "tsv_gene_differs_from_nearest_gtf_gene", False
    return True, "unclassified_conflict", False


def annotate_counts_df(counts_df: pd.DataFrame, args) -> pd.DataFrame:
    gtf_index = load_gtf_genes(args.gtf, args.gtf_feature, args.gene_name_attr, args.gene_id_attr)
    pos_df = counts_df[["contig", "position"]].drop_duplicates().copy()
    ann_rows = []
    for r in pos_df.itertuples(index=False):
        ann = annotate_pos(r.contig, int(r.position), gtf_index, args.slop, args.max_nearest_dist)
        ann_rows.append({"contig": r.contig, "position": int(r.position), **ann})
    ann_df = pd.DataFrame(ann_rows)
    out = counts_df.merge(ann_df, on=["contig", "position"], how="left")
    vals = [conflict_reason(g, a) for g, a in zip(out["gene"], out[["gtf_gene", "gtf_match_type"]].to_dict("records"))]
    out["gtf_gene_conflict"] = [v[0] for v in vals]
    out["gtf_conflict_reason"] = [v[1] for v in vals]
    out["tsv_gene_in_gtf_list"] = [v[2] for v in vals]
    return out


def sequential_filter_counts_df(annot_df: pd.DataFrame, args):
    if not (args.drop_conflicts or args.drop_multi_tsv_and_gtf):
        return annot_df, annot_df.iloc[0:0].copy()

    work = annot_df.copy()
    dropped_parts = []

    if args.drop_conflicts:
        mask = work["gtf_gene_conflict"].astype(bool)
        d = work.loc[mask].copy()
        d["drop_reason"] = "gtf_gene_conflict"
        dropped_parts.append(d)
        work = work.loc[~mask].copy()

    if args.drop_multi_tsv_and_gtf:
        pos_cols = ["contig", "position", "refAllele", "altAllele"]
        tmp = work.assign(gene_clean=work["gene"].map(norm_gene))
        tsv_counts = (
            tmp[tmp["gene_clean"] != ""]
            .groupby(pos_cols)["gene_clean"]
            .nunique()
            .reset_index(name="n_tsv_genes")
        )
        gtf_counts = work.groupby(pos_cols, as_index=False).agg(n_gtf_genes=("gtf_gene_n", "max"))
        flag = tsv_counts.merge(gtf_counts, on=pos_cols, how="outer").fillna(0)
        bad_pos = flag[(flag["n_tsv_genes"] >= 2) & (flag["n_gtf_genes"] >= 2)][pos_cols]
        if not bad_pos.empty:
            work2 = work.merge(bad_pos.assign(_drop_multi=True), on=pos_cols, how="left")
            mask = work2["_drop_multi"].eq(True)
            d = work2.loc[mask, work.columns].copy()
            d["drop_reason"] = "multi_tsv_gene_and_multi_gtf_gene"
            dropped_parts.append(d)
            work = work2.loc[~mask, work.columns].copy()

    dropped = pd.concat(dropped_parts, ignore_index=True) if dropped_parts else annot_df.iloc[0:0].copy()
    return work, dropped


def write_summary(summary_path: str, before_df: pd.DataFrame, final_df: pd.DataFrame, dropped_df: pd.DataFrame, args):
    pos_cols = ["contig", "position", "refAllele", "altAllele"]
    rows = [
        {"metric": "rows_total_before_filter", "value": len(before_df)},
        {"metric": "rows_written", "value": len(final_df)},
        {"metric": "rows_removed_total", "value": len(dropped_df)},
        {"metric": "drop_conflicts", "value": bool(args.drop_conflicts)},
        {"metric": "drop_multi_tsv_and_gtf", "value": bool(args.drop_multi_tsv_and_gtf)},
        {"metric": "unique_snp_positions_before_filter", "value": before_df[pos_cols].drop_duplicates().shape[0]},
        {"metric": "unique_snp_positions_written", "value": final_df[pos_cols].drop_duplicates().shape[0]},
        {"metric": "unique_snp_positions_removed", "value": dropped_df[pos_cols].drop_duplicates().shape[0] if not dropped_df.empty else 0},
    ]
    if "gtf_gene_conflict" in before_df.columns:
        rows.append({"metric": "rows_with_gtf_conflict_before_filter", "value": int(before_df["gtf_gene_conflict"].sum())})
    if "gtf_multi_gene" in before_df.columns:
        rows.append({"metric": "rows_with_multi_gene_overlap_or_nearby", "value": int(before_df["gtf_multi_gene"].fillna(False).sum())})
    if "gtf_conflict_reason" in before_df.columns:
        for k, v in before_df["gtf_conflict_reason"].value_counts(dropna=False).sort_index().items():
            rows.append({"metric": f"reason__{k}", "value": int(v)})
    if "drop_reason" in dropped_df.columns and not dropped_df.empty:
        for k, v in dropped_df["drop_reason"].value_counts(dropna=False).sort_index().items():
            rows.append({"metric": f"removed__{k}", "value": int(v)})
    pd.DataFrame(rows).to_csv(summary_path, sep="\t", index=False)


# --------------------- Core counting loop ---------------------

def load_sites_for_region(vcf, region: str, het_only: bool, sample_name: Optional[str]):
    site_order = []
    sites_by_pos0 = defaultdict(list)

    for rec in vcf.fetch(region=region):
        if not is_snp(rec):
            continue
        if het_only and not genotype_is_het_01(rec, sample_name):
            continue
        chrom = rec.contig
        pos1 = int(rec.pos)
        ref = rec.ref.upper()
        alt = rec.alts[0].upper()
        idx = len(site_order)
        site_order.append((chrom, pos1, ref, alt))
        sites_by_pos0[pos1 - 1].append(idx)

    return site_order, sites_by_pos0


def write_counts_for_site(out, site, counts_for_site) -> int:
    chrom, pos1, ref, alt = site
    rows = 0
    for (cell, gene), (rc, ac, oc) in counts_for_site.items():
        total = rc + ac + oc
        out.write(f"{cell}\t{chrom}\t{pos1}\t{gene}\t{ref}\t{alt}\t{rc}\t{ac}\t{total}\t{oc}\n")
        rows += 1
    return rows


def alignment_passes_filters(aln, args_dict) -> bool:
    if aln.is_unmapped:
        return False
    if not args_dict["allow_secondary"] and (aln.is_secondary or aln.is_supplementary):
        return False
    if not args_dict["count_duplicates"] and aln.is_duplicate:
        return False
    if aln.mapping_quality < int(args_dict["min_mapq"]):
        return False
    return True


def process_region_to_temp(bam_path: str, vcf_path: str, region: str, args_dict: dict, region_index: int):
    threads = int(args_dict["threads"])
    debug = bool(args_dict["debug"])
    debug_examples = int(args_dict["debug_examples"])

    bam = None
    vcf = None
    out = None
    ev_out = None
    tmp_path = None
    ev_tmp_path = None

    try:
        bam = pysam.AlignmentFile(bam_path, "rb", threads=threads)
        vcf = pysam.VariantFile(vcf_path, threads=threads)

        try:
            bam.check_index()
        except Exception as e:
            raise RuntimeError(f"[error] BAM index missing/unreadable for {bam_path}: {e}")

        het_only = bool(args_dict["het_only"])
        sample_name = get_sample_for_het_only(vcf, args_dict["sample"]) if het_only else None

        try:
            site_order, sites_by_pos0 = load_sites_for_region(vcf, region, het_only, sample_name)
        except Exception as e:
            raise RuntimeError(f"[error] VCF fetch failed for region {region}: {e}")

        fd, tmp_path = tempfile.mkstemp(prefix=f"counts_r{region_index:06d}_", suffix=".tsv")
        out = os.fdopen(fd, "w")
        out.write("\t".join(HEADER_COLS) + "\n")

        wants_evidence = bool(args_dict.get("read_evidence_out") or args_dict.get("dropped_read_ids_out") or args_dict.get("unique_dropped_read_ids_out") or args_dict.get("dropped_read_summary_out"))
        if wants_evidence:
            ev_fd, ev_tmp_path = tempfile.mkstemp(prefix=f"evidence_r{region_index:06d}_", suffix=".tsv")
            ev_out = os.fdopen(ev_fd, "w")
            ev_out.write("\t".join(EVIDENCE_COLS) + "\n")

        if len(site_order) == 0:
            if debug:
                print(f"[{region}] no SNPs to process", flush=True)
            return region_index, region, tmp_path, ev_tmp_path

        chrom, start0, end0 = region_to_fetch_args(region)
        flag_filter = build_flag_filter(args_dict["allow_secondary"], args_dict["count_duplicates"])
        cb_tag = args_dict["cb_tag"]
        skip_indels = bool(args_dict["skip_indels"])
        min_bq = int(args_dict["min_bq"])
        max_depth = int(args_dict["max_depth"])
        ignore_overlaps = bool(args_dict["ignore_overlaps"])

        counts_by_site = [defaultdict(lambda: [0, 0, 0]) for _ in site_order]
        examples_shown = 0
        pileup_cols_seen = 0

        try:
            col_iter = bam.pileup(
                chrom,
                start0,
                end0,
                truncate=True,
                stepper="samtools",
                max_depth=max_depth,
                flag_filter=flag_filter,
                ignore_overlaps=ignore_overlaps,
                min_base_quality=min_bq,
            )

            for col in col_iter:
                site_indexes = sites_by_pos0.get(col.reference_pos)
                if not site_indexes:
                    continue
                pileup_cols_seen += 1

                for pr in col.pileups:
                    if pr.is_refskip:
                        continue

                    aln = pr.alignment
                    if not alignment_passes_filters(aln, args_dict):
                        continue
                    if not aln.has_tag(cb_tag):
                        continue
                    cb = aln.get_tag(cb_tag)
                    if not cb:
                        continue
                    gene = aln.get_tag("GN") if aln.has_tag("GN") else "."

                    allele_class = None
                    debug_base = None

                    if pr.is_del:
                        if skip_indels:
                            continue
                        allele_class = 2
                        debug_base = "DEL"
                    elif pr.indel != 0:
                        if skip_indels:
                            continue
                        allele_class = 2
                        debug_base = f"INDEL({pr.indel})"
                    else:
                        qpos = pr.query_position
                        if qpos is None:
                            continue
                        try:
                            base = aln.query_sequence[qpos].upper()
                        except Exception:
                            continue
                        debug_base = base

                    for site_idx in site_indexes:
                        site = site_order[site_idx]
                        ref = site[2]
                        alt = site[3]

                        if allele_class is None:
                            if debug_base == ref:
                                allele_class_i = 0
                            elif debug_base == alt:
                                allele_class_i = 1
                            else:
                                allele_class_i = 2
                        else:
                            allele_class_i = allele_class

                        counts = counts_by_site[site_idx][(cb, gene)]
                        counts[allele_class_i] += 1

                        if ev_out is not None:
                            allele_label = "ref" if allele_class_i == 0 else ("alt" if allele_class_i == 1 else "other")
                            ev_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                aln.query_name, cb, site[0], site[1], gene, site[2], site[3], allele_label
                            ))

                    if debug and examples_shown < debug_examples:
                        print(f"[{region}] read={aln.query_name} CB={cb} GN={gene} base={debug_base}", flush=True)
                        examples_shown += 1

        except Exception as e:
            raise RuntimeError(f"[error] BAM pileup failed for region {region}: {e}")

        rows_written = 0
        for site_idx, site in enumerate(site_order):
            rows_written += write_counts_for_site(out, site, counts_by_site[site_idx])

        if debug:
            print(
                f"[{region}] done: SNPs={len(site_order)} pileup_site_columns={pileup_cols_seen} rows={rows_written}",
                flush=True,
            )

        return region_index, region, tmp_path, ev_tmp_path

    finally:
        if out is not None:
            out.close()
        if ev_out is not None:
            ev_out.close()
        if bam is not None:
            bam.close()
        if vcf is not None:
            vcf.close()


# ---------------------- Orchestration ----------------------

def merge_temp_files(results, out_path: str):
    results = sorted(results, key=lambda x: x[0])
    with open(out_path, "w") as w:
        w.write("\t".join(HEADER_COLS) + "\n")
        for item in results:
            tmp = item[2]
            with open(tmp) as f:
                next(f)
                for line in f:
                    w.write(line)


def merge_evidence_files(results, out_path: str):
    results = sorted(results, key=lambda x: x[0])
    with open(out_path, "w") as w:
        w.write("\t".join(EVIDENCE_COLS) + "\n")
        for item in results:
            ev_tmp = item[3] if len(item) > 3 else None
            if not ev_tmp:
                continue
            with open(ev_tmp) as f:
                next(f)
                for line in f:
                    w.write(line)


def cleanup_temp_files(results):
    for item in results:
        if not item:
            continue
        for tmp in item[2:]:
            try:
                if tmp and os.path.exists(tmp):
                    os.remove(tmp)
            except Exception:
                pass


def write_dropped_read_outputs(dropped_evidence: pd.DataFrame, args):
    if args.dropped_read_ids_out:
        dropped_evidence.to_csv(args.dropped_read_ids_out, sep="\t", index=False)
        print(f"[main] wrote dropped read x SNP evidence: {args.dropped_read_ids_out}", file=sys.stderr)

    if args.unique_dropped_read_ids_out:
        dropped_evidence[["read_id"]].drop_duplicates().to_csv(
            args.unique_dropped_read_ids_out,
            sep="\t",
            index=False,
            header=False,
        )
        print(f"[main] wrote unique dropped read IDs: {args.unique_dropped_read_ids_out}", file=sys.stderr)

    if args.dropped_read_summary_out:
        if dropped_evidence.empty:
            summary = pd.DataFrame(columns=[
                "read_id", "n_dropped_snp_rows", "n_dropped_snps", "cells",
                "genes", "drop_reasons", "allele_classes",
            ])
        else:
            summary = (
                dropped_evidence.groupby("read_id", as_index=False)
                .agg(
                    n_dropped_snp_rows=("position", "size"),
                    n_dropped_snps=("position", "nunique"),
                    cells=("cell_barcode", lambda s: ";".join(sorted(set(map(str, s))))),
                    genes=("gene", lambda s: ";".join(sorted(set(map(str, s))))),
                    drop_reasons=("drop_reason", lambda s: ";".join(sorted(set(map(str, s))))),
                    allele_classes=("allele_class", lambda s: ";".join(sorted(set(map(str, s))))),
                )
                .sort_values(["n_dropped_snps", "n_dropped_snp_rows", "read_id"], ascending=[False, False, True])
            )
        summary.to_csv(args.dropped_read_summary_out, sep="\t", index=False)
        print(f"[main] wrote dropped read summary: {args.dropped_read_summary_out}", file=sys.stderr)


def main(argv=None):
    args = parse_args(argv)
    regions = build_regions(args)

    if args.debug:
        print(f"[main] processing {len(regions)} region(s)", flush=True)

    results = []
    evidence_path = None
    raw_counts_path = None

    try:
        if args.procs <= 1:
            for i, region in enumerate(regions):
                results.append(process_region_to_temp(args.bam, args.vcf, region, vars(args), i))
        else:
            with ProcessPoolExecutor(max_workers=args.procs) as ex:
                futs = {
                    ex.submit(process_region_to_temp, args.bam, args.vcf, region, vars(args), i): (i, region)
                    for i, region in enumerate(regions)
                }
                for fut in as_completed(futs):
                    i, region = futs[fut]
                    try:
                        results.append(fut.result())
                    except Exception as e:
                        raise RuntimeError(f"Worker failed for region {region} (index {i}): {e}") from e

        needs_post = args.gtf is not None
        raw_counts_path = args.out if not needs_post else tempfile.mktemp(prefix="raw_counts_", suffix=".tsv")
        merge_temp_files(results, raw_counts_path)

        wants_any_evidence = bool(args.read_evidence_out or args.dropped_read_ids_out or args.unique_dropped_read_ids_out or args.dropped_read_summary_out)
        if wants_any_evidence:
            if args.read_evidence_out:
                evidence_path = args.read_evidence_out
            else:
                evidence_path = tempfile.mktemp(prefix="read_evidence_", suffix=".tsv")
            merge_evidence_files(results, evidence_path)

        if needs_post:
            counts_df = pd.read_csv(raw_counts_path, sep="\t")
            annot_df = annotate_counts_df(counts_df, args)
            final_df, dropped_df = sequential_filter_counts_df(annot_df, args)
            final_df.to_csv(args.out, sep="\t", index=False)

            summary_path = args.summary_out or args.out + ".summary.tsv"
            write_summary(summary_path, annot_df, final_df, dropped_df, args)
            print(f"[main] wrote summary: {summary_path}", file=sys.stderr)

            wants_dropped_outputs = bool(args.dropped_read_ids_out or args.unique_dropped_read_ids_out or args.dropped_read_summary_out)
            if wants_dropped_outputs:
                if dropped_df.empty:
                    empty = pd.DataFrame(columns=EVIDENCE_COLS + ["drop_reason"])
                    write_dropped_read_outputs(empty, args)
                else:
                    if evidence_path is None:
                        raise RuntimeError("Internal error: dropped read outputs requested but no evidence file was created")
                    pos_cols = ["cell_barcode", "contig", "position", "gene", "refAllele", "altAllele"]
                    dropped_keys = dropped_df[pos_cols + ["drop_reason"]].drop_duplicates()
                    ev = pd.read_csv(evidence_path, sep="\t")
                    ev2 = ev.merge(dropped_keys, on=pos_cols, how="inner")
                    write_dropped_read_outputs(ev2, args)

            if raw_counts_path != args.out and os.path.exists(raw_counts_path):
                os.remove(raw_counts_path)

        if evidence_path and evidence_path != args.read_evidence_out and os.path.exists(evidence_path):
            os.remove(evidence_path)

    finally:
        cleanup_temp_files(results)


if __name__ == "__main__":
    raise SystemExit(main())
