#!/usr/bin/env python3
"""
Haplotype phasing with:
  - Stage-1 connectivity check (save LCC vs non-LCC SNPs)
  - Expression-aware weighted graph partitioning + removal of same-partition SNPs
  - Stage-2 fill-in via per-cell active-X concordance
  - Validation metrics and outputs

CLI
---
phase_X_debug.py <filepath> <outdir>
                 [--min_reads 20] [--lower_cutoff 0.10] [--upper_cutoff None]
                 [--balance_tolerance IGNORED] [--repair_delta IGNORED]
                 [--discordant_tsv per_snp_concordance.tsv]

NB: --balance_tolerance and --repair_delta exist only so old commands don't crash;
    they are completely ignored in this version.
"""

import sys
import os
import argparse
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix
import numpy as np
import igraph as ig
from scipy.sparse import triu
import copy
import random
from typing import Dict, Tuple, List, Optional
import subprocess
# -----------------------
# CLI & setup
# -----------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Haplotype phasing with expression-aware weighted graph partitioning; repair/balancing disabled."
    )
    p.add_argument("filepath", help="Input TSV/CSV of allele counts")
    p.add_argument("outdir", help="Output directory")
    p.add_argument("--min_reads", type=int, default=20, help="Minimum total reads per SNP (default: 20)")
    p.add_argument("--lower_cutoff", type=float, default=0.10, help="Lower ref-fraction bound (default: 0.10)")
    p.add_argument(
        "--upper_cutoff",
        type=float,
        default=None,
        help="Upper ref-fraction bound (default: 1 - lower_cutoff)",
    )
    # kept for backward compatibility but ignored
    p.add_argument(
        "--balance_tolerance",
        type=int,
        default=None,
        help="IGNORED in this version (no partition balancing performed).",
    )
    p.add_argument(
        "--repair_delta",
        type=float,
        default=None,
        help="IGNORED in this version (no repair performed).",
    )
    p.add_argument(
        "--discordant_tsv",
        default=None,
        help="Optional per_snp_concordance.tsv to flag discordant SNPs and compute debug features",
    )
    p.add_argument(
        "--partition_mode",
        choices=["unweighted", "weighted", "expression"],
        default="expression",
        help=(
            "Graph partitioning mode. "
            "unweighted = old binary topology only; "
            "weighted = binary/scaled co-occurrence edge weights; "
            "expression = expression-aware cosine edge weights. "
            "Default: expression."
        ),
    )
    p.add_argument(
        "--edge_weight_mode",
        choices=["binary", "cosine", "overlap"],
        default=None,
        help=(
            "Advanced override for edge weights: binary, cosine, or overlap. "
            "If omitted, this is inferred from --partition_mode. "
            "Use --partition_mode for normal runs."
        ),
    )
    p.add_argument(
        "--vcf-out",
        default=None,
        help="Optional phased VCF output path. Default: <outdir>/haplotypes_min{min_reads}_{lc_tag}.vcf",
    )
    p.add_argument(
        "--sample-name",
        default="sample",
        help="Sample name to use in phased VCF. Default: sample",
    )
    p.add_argument(
        "--stage2-mode",
        choices=["cell", "expression"],
        default="cell",
        help=(
            "How to score Stage-2 fill-in SNPs. "
            "cell = original equal-cell vote based on ref>alt/alt>ref; "
            "expression = read-depth-weighted vote using allele-count margins. "
            "Default: cell."
        ),
    )
    p.add_argument(
        "--stage2-min-evidence",
        type=float,
        default=1.0,
        help=(
            "Minimum Stage-2 evidence required before assigning a SNP. "
            "For --stage2-mode cell this is informative cell count; "
            "for expression this is total informative allele reads. Default: 1."
        ),
    )
    p.add_argument(
        "--stage2-tie-action",
        choices=["random", "skip"],
        default="random",
        help="What to do when Stage-2 score is exactly 0.5. Default: random, matching old behaviour.",
    )
    return p.parse_args(argv)

def fmt_lc(x: float) -> str:
    return f"lc{float(x):.2f}"

# -----------------------
# I/O helpers
# -----------------------
def load_and_standardize(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.read_csv(path)

    cols_lower = {c.lower(): c for c in df.columns}
    new_schema = {
        "cell_barcode",
        "contig",
        "position",
        "gene",
        "refallele",
        "altallele",
        "refcount",
        "altcount",
        "totalcount",
        "othercount",
    }

    if new_schema.issubset(set(k.lower() for k in df.columns)):
        df = df.rename(
            columns={
                cols_lower["cell_barcode"]: "cell_barcode",
                cols_lower["position"]: "POS",
                cols_lower["refallele"]: "REF",
                cols_lower["altallele"]: "ALT",
                cols_lower["refcount"]: "REFcount",
                cols_lower["altcount"]: "ALTcount",
                cols_lower["gene"]: "Gene",
            }
        )
        if "contig" in cols_lower:
            df = df.rename(columns={cols_lower["contig"]: "CONTIG"})
        if "totalcount" in cols_lower:
            df = df.rename(columns={cols_lower["totalcount"]: "TOTALcount"})
        if "othercount" in cols_lower:
            df = df.rename(columns={cols_lower["othercount"]: "OTHcount"})
        if "SNP_ID" not in df.columns:
            df["SNP_ID"] = df.apply(
                lambda r: f"{str(r.get('CONTIG', 'chrX'))}:{r['POS']}:{r['REF']}:{r['ALT']}", axis=1
            )
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
        df = df.dropna(subset=["POS"])
        df["POS"] = df["POS"].astype(int)
        for col in ["REFcount", "ALTcount"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    else:
        if "POS" in df.columns:
            df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
            df = df.dropna(subset=["POS"])
            df["POS"] = df["POS"].astype(int)
        for col in ["REFcount", "ALTcount"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    if "Gene" in df.columns:
        def _blank(x):
            if pd.isna(x):
                return True
            s = str(x).strip()
            return s == "" or s.upper() == "NA" or s == "-"

        before = len(df)
        df = df[~df["Gene"].apply(_blank)].copy()
        after = len(df)
        print(f"Filtered blank Gene rows: {before - after} removed, {after} remain.")

    return df

# -----------------------
# Duplicate forensics
# -----------------------
def classify_dup_group(group):
    contigs = set(group["CONTIG"]) if "CONTIG" in group.columns else set()
    genes = set(group["Gene"]) if "Gene" in group.columns else set()
    alleles = (
        set(zip(group["REF"], group["ALT"]))
        if "REF" in group.columns and "ALT" in group.columns
        else set()
    )
    exact_dup = group.duplicated(keep=False).any()
    if len(alleles) > 1:
        return "multi-allelic"
    if len(genes) > 1:
        return "different-gene"
    if len(contigs) > 1:
        return "different-contig"
    if exact_dup:
        return "exact-duplicate-row"
    return "same-allele_repeated"

def debug_duplicates(df, outdir, examples_per_cause=10):
    if not {"cell_barcode", "POS"}.issubset(df.columns):
        print("Duplicate debug: missing required columns; skipping.")
        return
    dup_mask = df.duplicated(subset=["cell_barcode", "POS"], keep=False)
    dups = df[dup_mask].copy()
    if dups.empty:
        print("No duplicate (cell_barcode, POS) keys detected")
        return
    dups = dups.sort_values(["cell_barcode", "POS"])
    summaries = []
    for (cb, pos), grp in dups.groupby(["cell_barcode", "POS"], sort=False):
        cause = classify_dup_group(grp)
        summaries.append(
            {
                "cell_barcode": cb,
                "POS": pos,
                "n_rows": len(grp),
                "contigs_n": len(set(grp["CONTIG"])) if "CONTIG" in grp.columns else 0,
                "genes_n": len(set(grp["Gene"])) if "Gene" in grp.columns else 0,
                "allele_pairs_n": len(set(zip(grp["REF"], grp["ALT"])))
                if {"REF", "ALT"}.issubset(grp.columns)
                else 0,
                "cause": cause,
            }
        )
    summary_df = pd.DataFrame(summaries)
    cause_counts = summary_df["cause"].value_counts().sort_values(ascending=False)
    print("\nDuplicate (cell_barcode, POS) groups by cause:")
    for cause, n in cause_counts.items():
        print(f"  {cause:22s} : {n}")
    dups.to_csv(os.path.join(outdir, "duplicates_raw.tsv"), sep="\t", index=False)
    summary_df.to_csv(os.path.join(outdir, "duplicates_summary.tsv"), sep="\t", index=False)
    example_rows = []
    for cause in cause_counts.index:
        keys = (
            summary_df[summary_df["cause"] == cause]
            .head(examples_per_cause)[["cell_barcode", "POS"]]
            .to_records(index=False)
        )
        for cb, pos in keys:
            chunk = dups[(dups["cell_barcode"] == cb) & (dups["POS"] == pos)].copy()
            chunk["example_cause"] = cause
            example_rows.append(chunk)
    if example_rows:
        ex_df = pd.concat(example_rows, ignore_index=True)
        ex_df.to_csv(os.path.join(outdir, "duplicates_examples.tsv"), sep="\t", index=False)
    print(
        f"\nWrote duplicate diagnostics to:\n"
        f"  - {os.path.join(outdir, 'duplicates_raw.tsv')}\n"
        f"  - {os.path.join(outdir, 'duplicates_summary.tsv')}\n"
        f"  - {os.path.join(outdir, 'duplicates_examples.tsv')}"
    )

# -----------------------
# Matrix construction
# -----------------------
def load_ase(df):
    print("Formatting dataframe...")
    df = df.copy()
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df = df.dropna(subset=["POS"])
    df["POS"] = df["POS"].astype(int)
    for col in ["REFcount", "ALTcount"]:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    dup_pairs = df.duplicated(subset=["POS", "cell_barcode"], keep=False).sum()
    if dup_pairs:
        print(f"Note: aggregating {dup_pairs} duplicate rows by (POS, cell_barcode) via sum.")
    pos_ref_df = pd.pivot_table(
        df.assign(POS_ref=df["POS"].astype(str) + "_ref"),
        index="POS_ref",
        columns="cell_barcode",
        values="REFcount",
        aggfunc="sum",
        fill_value=0,
        observed=False,
    )
    pos_alt_df = pd.pivot_table(
        df.assign(POS_alt=df["POS"].astype(str) + "_alt"),
        index="POS_alt",
        columns="cell_barcode",
        values="ALTcount",
        aggfunc="sum",
        fill_value=0,
        observed=False,
    )
    formatted_df = pd.concat([pos_ref_df, pos_alt_df], axis=0).sort_index()
    formatted_df = formatted_df.apply(pd.to_numeric, errors="coerce").fillna(0)
    return formatted_df

# -----------------------
# Filter helpers
# -----------------------
def add_gene_column(df, pos_gene_dict):
    df = df.copy()
    df.loc[:, "Gene"] = df.index.str.split("_").str[0].astype(int).map(pos_gene_dict).astype("object")
    return df


def build_pos_gene_dict(df: pd.DataFrame) -> Dict[int, str]:
    """Choose one stable gene label per SNP position using read-count support.

    The input SNP count TSV may contain multiple rows for the same POS with
    different Gene values. A plain df.set_index("POS")["Gene"].to_dict() is
    order-dependent in that case. This function chooses the non-blank gene with
    the highest total REF+ALT support for each POS. Positions with no non-blank
    gene are omitted and later map to NA.
    """
    if "Gene" not in df.columns:
        return {}
    tmp = df.copy()
    tmp["Gene"] = tmp["Gene"].astype(str).str.strip()
    blank = {"", "-", ".", "NA", "NAN", "None", "nan"}
    tmp = tmp[~tmp["Gene"].isin(blank)].copy()
    if tmp.empty:
        return {}
    for col in ["REFcount", "ALTcount"]:
        tmp[col] = pd.to_numeric(tmp[col], errors="coerce").fillna(0)
    support = (
        tmp.groupby(["POS", "Gene"], as_index=False)[["REFcount", "ALTcount"]]
        .sum()
    )
    support["support"] = support["REFcount"] + support["ALTcount"]
    best = (
        support.sort_values(["POS", "support", "Gene"], ascending=[True, False, True])
        .drop_duplicates("POS", keep="first")
    )
    return best.set_index("POS")["Gene"].to_dict()

def filter_ase(formatted_df, pos_gene_dict, min_reads=5, lower_cutoff=0.01, upper_cutoff=0.99):
    split_index = formatted_df.index.str.split("_", expand=True)
    formatted_df.index = pd.MultiIndex.from_arrays(
        [split_index.get_level_values(0), split_index.get_level_values(1)], names=["ID", "type"]
    )
    summed_df = formatted_df.sum(axis=1).unstack("type")[["ref", "alt"]]
    summed_df["total"] = summed_df["ref"] + summed_df["alt"]
    summed_df["ref_fraction"] = summed_df["ref"] / summed_df["total"]
    filt = (summed_df["total"] >= min_reads) & summed_df["ref_fraction"].between(
        lower_cutoff, upper_cutoff, inclusive="both"
    )
    filtered_df = summed_df[filt]
    valid_ids = filtered_df.index.tolist()
    mask = formatted_df.index.get_level_values("ID").isin(valid_ids)
    filtered_formatted_df = formatted_df[mask]
    filtered_formatted_df.index = filtered_formatted_df.index.map(lambda x: f"{x[0]}_{x[1]}")
    num_snps = len(filtered_formatted_df.index.to_list()) // 2
    all_genes = set(add_gene_column(filtered_formatted_df, pos_gene_dict).Gene)
    print(f"\t{num_snps} SNP pairs identified in {len(all_genes)} genes.")
    original_filtered_formatted_df = copy.deepcopy(filtered_formatted_df)
    return filtered_formatted_df, all_genes, original_filtered_formatted_df

# -----------------------
# Graph helpers
# -----------------------
def compute_scaled_concordance(filtered_formatted_df, mode="cosine"):
    """
    Build an allele-by-allele concordance matrix.

    Rows are allele rows such as POS_ref / POS_alt; columns are cells.

    mode="binary" reproduces the old behaviour:
        weight(i,j) = n_cells(i and j) / min(n_cells(i), n_cells(j))

    mode="cosine" includes expression/read count information:
        weight(i,j) = dot(count_i, count_j) / (||count_i|| * ||count_j||)
    This keeps weights in [0,1] and downweights alleles that are merely observed
    together at very low or inconsistent expression.

    mode="overlap" also includes expression/read count information:
        weight(i,j) = sum_c min(count_i_c, count_j_c) / min(sum_c count_i_c, sum_c count_j_c)
    This is closest in spirit to the old containment-style normalization, but uses
    expression/read counts rather than only presence/absence. It can be slower for
    very dense matrices.
    """
    mode = str(mode).lower()
    m = csr_matrix(filtered_formatted_df.values.astype(np.float64))

    if mode == "binary":
        binary = m.astype(bool).astype(np.float64)
        C = binary.dot(binary.T)
        D = np.asarray(C.diagonal()).ravel()
        Cc = C.tocoo()
        rows, cols, data = Cc.row, Cc.col, Cc.data
        mins = np.minimum(D[rows], D[cols])
        valid = mins > 0
        scaled = np.zeros_like(data, dtype=np.float64)
        scaled[valid] = data[valid] / mins[valid]
        return coo_matrix((scaled, (rows, cols)), shape=C.shape).tocsr()

    if mode == "cosine":
        C = m.dot(m.T).tocoo()
        norms = np.sqrt(np.asarray(m.multiply(m).sum(axis=1)).ravel())
        denom = norms[C.row] * norms[C.col]
        valid = denom > 0
        scaled = np.zeros_like(C.data, dtype=np.float64)
        scaled[valid] = C.data[valid] / denom[valid]
        return coo_matrix((scaled, (C.row, C.col)), shape=C.shape).tocsr()

    if mode == "overlap":
        # Slower but interpretable expression-aware containment.
        dense = m.toarray()
        totals = dense.sum(axis=1)
        rows = []
        cols = []
        vals = []
        n = dense.shape[0]
        for i in range(n):
            xi = dense[i]
            if totals[i] <= 0:
                continue
            nz_i = xi > 0
            for j in range(i, n):
                if totals[j] <= 0:
                    continue
                xj = dense[j]
                if not np.any(nz_i & (xj > 0)):
                    continue
                denom = min(totals[i], totals[j])
                if denom <= 0:
                    continue
                val = np.minimum(xi, xj).sum() / denom
                if val > 0:
                    rows.extend([i, j] if i != j else [i])
                    cols.extend([j, i] if i != j else [j])
                    vals.extend([val, val] if i != j else [val])
        return coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    raise ValueError(f"Unknown edge_weight_mode: {mode}")

def construct_igraph(sparse_concordance_matrix):
    upper = triu(sparse_concordance_matrix, k=1, format="coo")

    # Do not create zero-weight edges. This matters because igraph treats an edge
    # with weight 0 as still topologically present unless removed.
    keep = upper.data > 0
    rows = upper.row[keep]
    cols = upper.col[keep]
    weights = upper.data[keep].astype(float)

    g = ig.Graph(n=sparse_concordance_matrix.shape[0])
    if len(rows):
        g.add_edges(list(zip(rows, cols)))
        g.es["weight"] = weights.tolist()
    else:
        g.es["weight"] = []

    # Remove the ref-alt edge for each SNP (if present).
    # These two alleles are mutually exclusive biological alternatives and should
    # not directly pull each other into the same community.
    num_snps = sparse_concordance_matrix.shape[0] // 2
    pairs = np.column_stack((2 * np.arange(num_snps), 2 * np.arange(num_snps) + 1))
    edge_ids = g.get_eids(pairs.tolist(), directed=False, error=False)
    existing = np.array(edge_ids)[np.array(edge_ids) != -1]
    if existing.size > 0:
        g.delete_edges(existing.tolist())
    return g

def split_graph(g, use_weights: bool = True):
    """2-way graph split, optionally using edge weights.

    use_weights=False reproduces the old unweighted topology-only partition.
    use_weights=True uses g.es["weight"] when available.
    """
    if g.vcount() == 0:
        return []
    if g.ecount() == 0:
        # Degenerate case: no information. Split deterministically by vertex index
        # so downstream code does not crash, but these results should be treated as
        # uninformative.
        return [i % 2 for i in range(g.vcount())]

    weights = g.es["weight"] if (use_weights and "weight" in g.es.attributes() and g.ecount()) else None

    try:
        part = g.community_leading_eigenvector(weights=weights, clusters=2)
    except TypeError:
        part = g.community_leading_eigenvector(weights=weights)

    mem = np.asarray(part.membership, dtype=int)
    labels = sorted(set(mem.tolist()))

    if len(labels) == 1:
        return [0 for _ in mem]

    if len(labels) == 2:
        label_map = {labels[0]: 0, labels[1]: 1}
        return [label_map[x] for x in mem]

    # Fallback for igraph versions that return >2 communities despite no clusters arg.
    # This preserves a 2-haplotype model without silently treating community 2/3/etc
    # as arbitrary extra partitions.
    sizes = pd.Series(mem).value_counts()
    largest = sizes.index[0]
    mem2 = np.where(mem == largest, 0, 1)
    print(f"	[WARN] leading_eigenvector returned {len(labels)} communities; collapsed to largest-vs-rest")
    return mem2.tolist()

# -----------------------
# Connectivity check
# -----------------------
def run_stage1_connectivity_check(df_std, pos_gene_dict, min_reads, lower_cutoff, upper_cutoff, outdir, lc_tag, edge_weight_mode="cosine"):
    """
    Builds filtered allele matrix, constructs the Stage-1 graph, checks connectivity,
    and writes LCC vs non-LCC SNP lists.
    """
    formatted_df = load_ase(df_std)
    filtered_formatted_df, _, _ = filter_ase(
        formatted_df, pos_gene_dict, min_reads=min_reads, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff
    )
    print("[CONNECTIVITY] Building Stage-1 graph for connectivity check...")
    sc = compute_scaled_concordance(filtered_formatted_df, mode=edge_weight_mode)
    g = construct_igraph(sc)

    if g.vcount() == 0:
        print("[CONNECTIVITY] Empty graph after filtering.")
        lcc_out = os.path.join(outdir, f"stage1_lcc_snps_min{min_reads}_{lc_tag}.tsv")
        non_out = os.path.join(outdir, f"stage1_non_lcc_snps_min{min_reads}_{lc_tag}.tsv")
        pd.DataFrame(columns=["POS", "status"]).to_csv(lcc_out, sep="\t", index=False)
        pd.DataFrame(columns=["POS", "status"]).to_csv(non_out, sep="\t", index=False)
        return

    comps = g.components(mode="WEAK")
    comp_id = np.zeros(g.vcount(), dtype=int)
    for idx, comp in enumerate(comps):
        comp_id[comp] = idx

    names = filtered_formatted_df.index.to_list()
    pos_ref_comp = {}
    pos_alt_comp = {}
    for i in range(0, len(names), 2):
        pos = int(names[i].split("_")[0])
        pos_ref_comp[pos] = comp_id[i]
        pos_alt_comp[pos] = comp_id[i + 1]

    comp_sizes_vertices = [len(c) for c in comps]
    lcc_id = int(np.argmax(comp_sizes_vertices))
    rows = []
    for pos in sorted(pos_ref_comp.keys()):
        cr = pos_ref_comp[pos]
        ca = pos_alt_comp[pos]
        if cr == lcc_id and ca == lcc_id:
            status = "lcc"
        elif cr != ca:
            status = "split"
        else:
            status = "non_lcc"
        rows.append({"POS": pos, "comp_ref": int(cr), "comp_alt": int(ca), "status": status})

    comp_df = pd.DataFrame(rows).sort_values("POS")
    if {"REF", "ALT"}.issubset(df_std.columns):
        ref_lookup = df_std.drop_duplicates("POS").set_index("POS")["REF"].to_dict()
        alt_lookup = df_std.drop_duplicates("POS").set_index("POS")["ALT"].to_dict()
        comp_df["REF"] = comp_df["POS"].map(ref_lookup)
        comp_df["ALT"] = comp_df["POS"].map(alt_lookup)
    if "Gene" in df_std.columns:
        gene_lookup = df_std.drop_duplicates("POS").set_index("POS")["Gene"].to_dict()
        comp_df["Gene"] = comp_df["POS"].map(gene_lookup)

    lcc_snps = comp_df[comp_df["status"] == "lcc"].copy()
    non_lcc_snps = comp_df[comp_df["status"] != "lcc"].copy()

    lcc_out = os.path.join(outdir, f"stage1_lcc_snps_min{min_reads}_{lc_tag}.tsv")
    non_out = os.path.join(outdir, f"stage1_non_lcc_snps_min{min_reads}_{lc_tag}.tsv")
    lcc_snps.to_csv(lcc_out, sep="\t", index=False)
    non_lcc_snps.to_csv(non_out, sep="\t", index=False)

    n_comp = len(comps)
    n_pairs = len(names) // 2
    print(
        f"[CONNECTIVITY] Components: {n_comp} | LCC pairs: {len(lcc_snps)} / {n_pairs} "
        f"({len(lcc_snps)/max(1,n_pairs):.1%}); non-LCC (incl split): {len(non_lcc_snps)}"
    )
    print(f"[CONNECTIVITY] Wrote:\n  - {lcc_out}\n  - {non_out}")

# -----------------------
# Partition helpers & iterative removal
# -----------------------
def _count_snps_per_partition(partition_array: List[int]) -> Tuple[int, int]:
    p0 = p1 = 0
    for i in range(len(partition_array) // 2):
        ref_part = partition_array[2 * i]
        if ref_part == 0:
            p0 += 1
        else:
            p1 += 1
    return p0, p1

def remove_duplicate_SNPs(filtered_formatted_df, partition_array):
    to_drop_idx = []
    removed_p0 = 0
    removed_p1 = 0
    for i in range(len(partition_array) // 2):
        r = partition_array[2 * i]
        a = partition_array[2 * i + 1]
        if r == a:
            to_drop_idx += [2 * i, 2 * i + 1]
            if r == 0:
                removed_p0 += 1
            else:
                removed_p1 += 1
    n_removed = len(to_drop_idx) // 2
    removable = np.array(filtered_formatted_df.index.to_list())[to_drop_idx]
    filtered_formatted_df = filtered_formatted_df[~filtered_formatted_df.index.isin(removable)]
    return filtered_formatted_df, n_removed, {"P0": removed_p0, "P1": removed_p1}

def iteratively_remove_SNPs(
    filtered_formatted_df,
    log_path: Optional[str] = None,
    edge_weight_mode: str = "cosine",
    use_partition_weights: bool = True,
):
    """
    Iteratively:
      - build graph
      - partition (leading eigenvector)
      - drop all SNPs with ref/alt in same partition
    until no removals.

    NO repair and NO artificial haplotype balancing.
    """
    print("Constructing graph...")
    n_removed = -1
    it = 0
    g = None
    part = None
    log_rows = []

    while n_removed != 0:
        sc = compute_scaled_concordance(filtered_formatted_df, mode=edge_weight_mode)
        g = construct_igraph(sc)
        part = split_graph(g, use_weights=use_partition_weights)

        p0_snps, p1_snps = _count_snps_per_partition(part)
        total_pairs_now = len(part) // 2

        # drop all same-partition pairs
        filtered_formatted_df, n_removed, removed_counts = remove_duplicate_SNPs(
            filtered_formatted_df, part
        )
        removed_p0 = removed_counts["P0"]
        removed_p1 = removed_counts["P1"]

        it += 1
        remaining = len(filtered_formatted_df.index.to_list()) // 2

        log_rows.append(
            {
                "iteration": it,
                "snps_total": total_pairs_now,
                "p0_snps": p0_snps,
                "p1_snps": p1_snps,
                "repaired_total": 0,
                "repaired_to_p0": 0,
                "repaired_to_p1": 0,
                "removed_total": n_removed,
                "removed_p0": removed_p0,
                "removed_p1": removed_p1,
                "remaining_total": remaining,
                "cross_pairs": total_pairs_now - n_removed,
                "removed_frac": (n_removed / total_pairs_now) if total_pairs_now else float("nan"),
            }
        )

        print(
            f"\tPartition before removal: P0={p0_snps}, P1={p1_snps} (pairs={total_pairs_now}); "
            f"removed={n_removed} (removed P0={removed_p0}, P1={removed_p1}); remaining={remaining}"
        )

        if n_removed == 0:
            break

    left = len(filtered_formatted_df.index.to_list()) // 2
    print(f"\tConverged after {it} iterations. {left} SNP pairs left.")

    if log_path is not None and log_rows:
        pd.DataFrame(log_rows).to_csv(log_path, sep="\t", index=False)
        print(f"\t[LOG] Wrote Stage-1 partition/removal log: {log_path}")

    return filtered_formatted_df, g, part

# -----------------------
# Stage-1 → haplotypes
# -----------------------
def generate_haplotype(
    df,
    pos_gene_dict,
    min_reads=5,
    lower_cutoff=0.0,
    upper_cutoff=1.0,
    log_path: Optional[str] = None,
    edge_weight_mode: str = "cosine",
    use_partition_weights: bool = True,
):
    formatted_df = load_ase(df)
    filtered_formatted_df, all_genes, original_filtered_formatted_df = filter_ase(
        formatted_df, pos_gene_dict, min_reads=min_reads, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff
    )

    filtered_formatted_df, g, partition_array = iteratively_remove_SNPs(
        filtered_formatted_df,
        log_path=log_path,
        edge_weight_mode=edge_weight_mode,
        use_partition_weights=use_partition_weights,
    )

    p0_final, p1_final = _count_snps_per_partition(partition_array)
    print(
        f"\tFinal partitions (no rebalancing): "
        f"P0={p0_final}, P1={p1_final}"
    )

    p1_idx = [i for i, x in enumerate(partition_array) if x == 0]
    p2_idx = [i for i, x in enumerate(partition_array) if x == 1]
    hap1 = np.array(filtered_formatted_df.index.to_list())[p1_idx]
    hap2 = np.array(filtered_formatted_df.index.to_list())[p2_idx]

    # No artificial phase balancing.
    # Each retained SNP should contribute exactly one allele to H1 and one allele to H2
    # if the graph split placed ref/alt in opposite communities. Same-community SNPs
    # have already been removed by iteratively_remove_SNPs().
    hap1 = hap1.tolist() if isinstance(hap1, np.ndarray) else list(hap1)
    hap2 = hap2.tolist() if isinstance(hap2, np.ndarray) else list(hap2)
    return filtered_formatted_df, hap1, hap2, g, partition_array, all_genes, original_filtered_formatted_df

# -----------------------
# Stage 2 and validation helpers
# -----------------------
def calculate_Xi(filepath, filtered_formatted_df, haplotype_1, haplotype_2, pos_gene_dict):
    # Use stable support-aware POS->Gene mapping built once in main.
    POS_Gene = pos_gene_dict
    h1c = filtered_formatted_df.loc[filtered_formatted_df.index.isin(haplotype_1)]
    h2c = filtered_formatted_df.loc[filtered_formatted_df.index.isin(haplotype_2)]

    def add_gene_column_inner(df):
        df = df.copy()
        df.loc[:, "Gene"] = df.index.str.split("_").str[0].astype(int).map(POS_Gene).astype("object")
        return df

    h1c = add_gene_column_inner(h1c)
    h2c = add_gene_column_inner(h2c)

    def collapse(df, how="max"):
        return df.groupby("Gene").max() if how == "max" else df.groupby("Gene").mean()

    c1 = collapse(h1c, "max")
    c2 = collapse(h2c, "max")
    s1 = c1.sum()
    s2 = c2.sum()
    sums = pd.DataFrame({"Haplotype_1_Sum": s1, "Haplotype_2_Sum": s2})
    sums["active_X"] = sums.apply(
        lambda r: "X1"
        if r["Haplotype_1_Sum"] > r["Haplotype_2_Sum"]
        else ("X2" if r["Haplotype_1_Sum"] < r["Haplotype_2_Sum"] else "X1/X2"),
        axis=1,
    )
    return h1c, h2c, sums, sums["active_X"].value_counts(), c1, c2

def phase_missing_snps(
    original_filtered_formatted_df, filtered_formatted_df, haplotype_1, haplotype_2, c1, c2, sums
):
    original_snps = original_filtered_formatted_df.index.to_list()
    discarded = [s for s in original_snps if s not in haplotype_1 and s not in haplotype_2]
    comp = {}
    for cell in original_filtered_formatted_df.columns:
        res = {}
        for i in range(len(discarded) // 2):
            s1 = discarded[2 * i]
            s2 = discarded[2 * i + 1]
            name = s1.split("_")[0]
            if (
                s1 in original_filtered_formatted_df.index
                and s2 in original_filtered_formatted_df.index
            ):
                v1 = original_filtered_formatted_df.at[s1, cell]
                v2 = original_filtered_formatted_df.at[s2, cell]
                res[name] = "X1" if v1 > v2 else ("X2" if v1 < v2 else "X1/X2")
            else:
                res[name] = "NA"
        comp[cell] = res
    return pd.DataFrame.from_dict(comp, orient="index")

def wilson_ci(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    if n <= 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z**2 / n
    center = p + z**2 / (2 * n)
    margin = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n)
    lo = (center - margin) / denom
    hi = (center + margin) / denom
    return (float(lo), float(hi))

def compute_partition_validation(
    g: ig.Graph, part: List[int], filtered_formatted_df: pd.DataFrame
) -> pd.DataFrame:
    if g is None or g.vcount() == 0:
        return (
            pd.DataFrame(
                columns=[
                    "POS",
                    "graph_deg_ref",
                    "graph_deg_alt",
                    "graph_deg_sum",
                    "graph_strength_ref",
                    "graph_strength_alt",
                    "graph_strength_sum",
                    "mean_edge_weight_ref",
                    "mean_edge_weight_alt",
                    "mean_edge_weight_avg",
                    "within_w_sum",
                    "cross_w_sum",
                    "partition_separation_w",
                    "within_deg",
                    "cross_deg",
                    "partition_separation_deg",
                ]
            )
            .set_index("POS")
        )

    names = filtered_formatted_df.index.to_list()
    pos2pair = {}
    for i in range(0, len(names), 2):
        p = int(names[i].split("_")[0])
        pos2pair[p] = (i, i + 1)

    deg = g.degree()
    try:
        strength = g.strength(weights=g.es["weight"])
    except Exception:
        strength = deg

    incident_edges = [g.incident(v) for v in range(g.vcount())]
    weights = g.es["weight"] if g.ecount() else []

    rows = []
    for pos, (i_ref, i_alt) in pos2pair.items():

        def allele_stats(i):
            degu = deg[i]
            strw = strength[i] if strength is not None else np.nan
            inc = incident_edges[i]
            if len(inc) == 0:
                return degu, strw, np.nan, 0.0, 0.0, 0, 0
            w = np.array([weights[e] for e in inc], dtype=float)
            mean_w = float(np.mean(w)) if len(w) else np.nan
            myp = part[i]
            nbrs = [
                g.es[e].tuple[1] if g.es[e].tuple[0] == i else g.es[e].tuple[0]
                for e in inc
            ]
            nbrs = np.array(nbrs, dtype=int)
            same = part if isinstance(part, np.ndarray) else np.array(part)
            same_mask = same[nbrs] == myp
            within_w = float(w[same_mask].sum()) if same_mask.any() else 0.0
            cross_w = float(w[~same_mask].sum()) if (~same_mask).any() else 0.0
            within_d = int(same_mask.sum())
            cross_d = int((~same_mask).sum())
            return degu, strw, mean_w, within_w, cross_w, within_d, cross_d

        dref, sref, mwref, winw_ref, crossw_ref, wind_ref, crossd_ref = allele_stats(i_ref)
        da, sa, mwa, winw_alt, crossw_alt, wind_alt, crossd_alt = allele_stats(i_alt)

        graph_deg_sum = int(dref + da)
        graph_strength_sum = float(
            (sref if sref == sref else 0) + (sa if sa == sa else 0)
        )
        mean_edge_weight_avg = np.nanmean([mwref, mwa])

        within_w_sum = winw_ref + winw_alt
        cross_w_sum = crossw_ref + crossw_alt
        sep_w = (
            (within_w_sum - cross_w_sum) / (within_w_sum + cross_w_sum)
            if (within_w_sum + cross_w_sum) > 0
            else np.nan
        )

        within_deg = wind_ref + wind_alt
        cross_deg = crossd_ref + crossd_alt
        sep_d = (
            (within_deg - cross_deg) / (within_deg + cross_deg)
            if (within_deg + cross_deg) > 0
            else np.nan
        )

        rows.append(
            {
                "POS": pos,
                "graph_deg_ref": int(dref),
                "graph_deg_alt": int(da),
                "graph_deg_sum": graph_deg_sum,
                "graph_strength_ref": float(sref) if sref == sref else np.nan,
                "graph_strength_alt": float(sa) if sa == sa else np.nan,
                "graph_strength_sum": graph_strength_sum,
                "mean_edge_weight_ref": float(mwref) if mwref == mwref else np.nan,
                "mean_edge_weight_alt": float(mwa) if mwa == mwa else np.nan,
                "mean_edge_weight_avg": float(mean_edge_weight_avg)
                if mean_edge_weight_avg == mean_edge_weight_avg
                else np.nan,
                "within_w_sum": float(within_w_sum),
                "cross_w_sum": float(cross_w_sum),
                "partition_separation_w": float(sep_w) if sep_w == sep_w else np.nan,
                "within_deg": int(within_deg),
                "cross_deg": int(cross_deg),
                "partition_separation_deg": float(sep_d) if sep_d == sep_d else np.nan,
            }
        )

    return pd.DataFrame(rows).set_index("POS")

def compute_xa_xi_reads(original_filtered_formatted_df, haplotype_1_all, haplotype_2_all, active_calls):
    if isinstance(active_calls, pd.DataFrame):
        active_calls = active_calls["active_X"]
    active_calls = active_calls.astype(str)
    active_calls = active_calls.reindex(original_filtered_formatted_df.columns)

    h1_set = set(haplotype_1_all)
    h2_set = set(haplotype_2_all)
    has_row = set(original_filtered_formatted_df.index)
    out = {}

    mask_x1 = active_calls.eq("X1").values
    mask_x2 = active_calls.eq("X2").values

    positions = sorted({int(x.split("_")[0]) for x in (h1_set | h2_set)})
    for pos in positions:
        rk = f"{pos}_ref"
        ak = f"{pos}_alt"
        if rk not in has_row or ak not in has_row:
            continue

        if rk in h1_set:
            h1_reads = original_filtered_formatted_df.loc[rk].values
            h2_reads = original_filtered_formatted_df.loc[ak].values
        elif ak in h1_set:
            h1_reads = original_filtered_formatted_df.loc[ak].values
            h2_reads = original_filtered_formatted_df.loc[rk].values
        else:
            continue

        xa = (h1_reads[mask_x1].sum() if mask_x1.any() else 0) + (
            h2_reads[mask_x2].sum() if mask_x2.any() else 0
        )
        xi = (h2_reads[mask_x1].sum() if mask_x1.any() else 0) + (
            h1_reads[mask_x2].sum() if mask_x2.any() else 0
        )

        out[pos] = {"Xa_reads": int(xa), "Xi_reads": int(xi)}

    return out

def compute_haplotype_consistency(
    original_filtered_formatted_df, haplotype_1_all, haplotype_2_all
) -> pd.DataFrame:
    h1_set = set(haplotype_1_all)
    h2_set = set(haplotype_2_all)
    has_row = set(original_filtered_formatted_df.index)

    rows = []
    positions = sorted({int(x.split("_")[0]) for x in (h1_set | h2_set)})
    for pos in positions:
        rk = f"{pos}_ref"
        ak = f"{pos}_alt"
        if rk not in has_row or ak not in has_row:
            continue

        if rk in h1_set:
            h1 = original_filtered_formatted_df.loc[rk].astype(float).values
            h2 = original_filtered_formatted_df.loc[ak].astype(float).values
        else:
            h1 = original_filtered_formatted_df.loc[ak].astype(float).values
            h2 = original_filtered_formatted_df.loc[rk].astype(float).values

        tot = h1 + h2
        mask = tot > 0
        if mask.sum() == 0:
            rows.append(
                {
                    "POS": pos,
                    "H1_wins_frac": np.nan,
                    "n_informative_cells_phase": 0,
                    "H1_wins_frac_ci_lo": np.nan,
                    "H1_wins_frac_ci_hi": np.nan,
                }
            )
            continue

        wins = h1[mask] > h2[mask]
        ties = h1[mask] == h2[mask]
        info_mask = ~ties
        n = int(info_mask.sum())
        if n == 0:
            rows.append(
                {
                    "POS": pos,
                    "H1_wins_frac": np.nan,
                    "n_informative_cells_phase": 0,
                    "H1_wins_frac_ci_lo": np.nan,
                    "H1_wins_frac_ci_hi": np.nan,
                }
            )
            continue
        k = int(wins[info_mask].sum())
        frac = k / n
        lo, hi = wilson_ci(k, n, z=1.96)
        rows.append(
            {
                "POS": pos,
                "H1_wins_frac": float(frac),
                "n_informative_cells_phase": n,
                "H1_wins_frac_ci_lo": float(lo),
                "H1_wins_frac_ci_hi": float(hi),
            }
        )

    return pd.DataFrame(rows).set_index("POS")

def _compute_pos_features(original_filtered_formatted_df, pos_list):
    features = []
    all_cells = list(original_filtered_formatted_df.columns)
    n_cells = len(all_cells)

    have = set(original_filtered_formatted_df.index)
    for pos in pos_list:
        ref_key = f"{pos}_ref"
        alt_key = f"{pos}_alt"
        if (ref_key not in have) or (alt_key not in have):
            features.append(
                {
                    "POS": int(pos),
                    "total_ref_sum": np.nan,
                    "total_alt_sum": np.nan,
                    "total_coverage": np.nan,
                    "n_cells_nonzero": np.nan,
                    "frac_cells_nonzero": np.nan,
                    "cell_ref_fraction_median": np.nan,
                    "cell_ref_fraction_iqr": np.nan,
                    "cell_ref_fraction_std": np.nan,
                }
            )
            continue

        r = original_filtered_formatted_df.loc[ref_key].astype(float).values
        a = original_filtered_formatted_df.loc[alt_key].astype(float).values
        tot = r + a

        total_ref_sum = float(r.sum())
        total_alt_sum = float(a.sum())
        total_cov = float(tot.sum())

        mask = tot > 0
        n_nonzero = int(mask.sum())
        frac_nonzero = n_nonzero / n_cells if n_cells else np.nan

        if n_nonzero > 0:
            ref_frac_cells = r[mask] / tot[mask]
            med = float(np.median(ref_frac_cells))
            q75 = float(np.percentile(ref_frac_cells, 75))
            q25 = float(np.percentile(ref_frac_cells, 25))
            iqr = q75 - q25
            std = float(np.std(ref_frac_cells))
        else:
            med = iqr = std = np.nan

        features.append(
            {
                "POS": int(pos),
                "total_ref_sum": total_ref_sum,
                "total_alt_sum": total_alt_sum,
                "total_coverage": total_cov,
                "n_cells_nonzero": n_nonzero,
                "frac_cells_nonzero": frac_nonzero,
                "cell_ref_fraction_median": med,
                "cell_ref_fraction_iqr": iqr,
                "cell_ref_fraction_std": std,
            }
        )

    return pd.DataFrame(features).set_index("POS")

def compute_pos_level_features(original_filtered_formatted_df, pos_list: List[int]) -> pd.DataFrame:
    return _compute_pos_features(original_filtered_formatted_df, pos_list)

def assemble_validation_matrix(
    filtered_formatted_df: pd.DataFrame,
    original_filtered_formatted_df: pd.DataFrame,
    g: Optional[ig.Graph],
    part: Optional[List[int]],
    haplotype_1_all: List[str],
    haplotype_2_all: List[str],
    active_calls: pd.Series,
    xa_xi_reads: Dict[int, Dict[str, int]],
) -> pd.DataFrame:
    val_graph = (
        compute_partition_validation(g, part, filtered_formatted_df)
        if (g is not None and part is not None)
        else pd.DataFrame()
    )
    val_cons = compute_haplotype_consistency(
        original_filtered_formatted_df, haplotype_1_all, haplotype_2_all
    )
    positions_all = sorted(
        {int(x.split("_")[0]) for x in (set(haplotype_1_all) | set(haplotype_2_all))}
    )
    val_cov = compute_pos_level_features(original_filtered_formatted_df, positions_all)

    xa_xi_df = pd.DataFrame.from_dict(xa_xi_reads, orient="index")
    xa_xi_df.index.name = "POS"
    xa_xi_df["xa_xi_margin"] = (xa_xi_df["Xa_reads"] - xa_xi_df["Xi_reads"]) / xa_xi_df[
        "Xa_reads"
    ].add(xa_xi_df["Xi_reads"]).replace({0: np.nan})

    parts = []
    if not val_graph.empty:
        parts.append(val_graph)
    parts += [val_cons, val_cov, xa_xi_df]
    val = parts[0].join(parts[1:], how="outer") if parts else pd.DataFrame(index=positions_all)

    def center01(x):  # [-1,1] -> [0,1]
        return 0.5 * (x + 1) if pd.notnull(x) else np.nan

    conf_components = []
    if "partition_separation_w" in val:
        conf_components.append(val["partition_separation_w"].map(center01))
    if "H1_wins_frac" in val:
        conf_components.append((val["H1_wins_frac"] - 0.5).abs() * 2)
    if "xa_xi_margin" in val:
        conf_components.append(val["xa_xi_margin"].abs())
    if conf_components:
        comp_mat = pd.concat(conf_components, axis=1)
        val["confidence_basic"] = comp_mat.mean(axis=1, skipna=True)

    return val

# -----------------------
# Minimal per-SNP haplotype table
# -----------------------
def build_min_haplotype_table(
    df_std,
    haplotype_1_all,
    haplotype_2_all,
    pos_gene_dict,
    stage_map,
    coverage_map,
):
    h1_set = set(haplotype_1_all)
    h2_set = set(haplotype_2_all)
    positions = sorted({int(s.split("_")[0]) for s in (h1_set | h2_set)})

    contig_lookup = (
        df_std.loc[:, ["POS", "CONTIG"]]
        .drop_duplicates("POS")
        .set_index("POS")["CONTIG"]
        .to_dict()
    ) if "CONTIG" in df_std.columns else {}

    ref_lookup = (
        df_std.loc[:, ["POS", "REF"]]
        .drop_duplicates("POS")
        .set_index("POS")["REF"]
        .to_dict()
    )
    alt_lookup = (
        df_std.loc[:, ["POS", "ALT"]]
        .drop_duplicates("POS")
        .set_index("POS")["ALT"]
        .to_dict()
    )

    rows = []
    for pos in positions:
        ref_base = ref_lookup.get(pos, None)
        alt_base = alt_lookup.get(pos, None)
        if ref_base is None or alt_base is None:
            continue

        ref_key = f"{pos}_ref"
        alt_key = f"{pos}_alt"

        if ref_key in h1_set:
            h1_base, h2_base = ref_base, alt_base
        elif alt_key in h1_set:
            h1_base, h2_base = alt_base, ref_base
        else:
            continue

        st = stage_map.get(pos, {})
        stage_val = float(st.get("stage")) if st.get("stage") is not None else float("nan")
        method = st.get("method", None)
        cov = coverage_map.get(pos, float("nan"))
        cov = int(cov) if pd.notnull(cov) else cov

        rows.append(
            {
                "CONTIG": contig_lookup.get(pos, "chrX"),
                "POS": int(pos),
                "REF": ref_base,
                "ALT": alt_base,
                "Gene": pos_gene_dict.get(pos, None),
                "H1": h1_base,
                "H2": h2_base,
                "Stage": stage_val,
                "Stage2_method": "" if method is None else method,
                "Total_coverage": cov,
            }
        )

    out_df = (
        pd.DataFrame(rows)
        .sort_values("POS", kind="mergesort")
        .loc[:, ["CONTIG", "POS", "REF", "ALT", "Gene", "H1", "H2", "Stage", "Stage2_method", "Total_coverage"]]
    )

    # Xist is transcribed from the inactive X, opposite to the usual
    # active-X-derived orientation used for the other genes. Flip its final
    # H1/H2 assignment only in the output table/VCF; do not alter the graph
    # partitioning, active-X calls, or Stage-2 calculations upstream.
    if not out_df.empty:
        xist_mask = out_df["Gene"].astype(str).str.lower().eq("xist")
        n_xist = int(xist_mask.sum())
        if n_xist:
            out_df.loc[xist_mask, ["H1", "H2"]] = out_df.loc[xist_mask, ["H2", "H1"]].to_numpy()
            print(f"[XIST] Flipped final H1/H2 assignment for {n_xist} Xist SNPs")

    return out_df


# -----------------------
# Phased VCF output
# -----------------------
def write_phased_vcf(min_df: pd.DataFrame, out_path: str, sample_name: str = "sample"):
    """Write a minimal phased VCF from the haplotype table.

    Encoding convention:
      - GT = 0|1 when H1 == REF and H2 == ALT
      - GT = 1|0 when H1 == ALT and H2 == REF

    This makes phase-set/sample haplotypes explicit while preserving the input
    REF/ALT alleles. Stage and gene info are stored in INFO.
    """
    if min_df.empty:
        raise ValueError("Cannot write VCF: haplotype table is empty")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    df = min_df.copy().sort_values("POS", kind="mergesort")

    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=phase_X_debug\n")
        out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Support-aware gene label for SNP">\n')
        out.write('##INFO=<ID=STAGE,Number=1,Type=String,Description="Phasing stage: 1 graph, 2 fill-in">\n')
        out.write('##INFO=<ID=STAGE2_METHOD,Number=1,Type=String,Description="Stage 2 method if applicable">\n')
        out.write('##INFO=<ID=TOTAL_COVERAGE,Number=1,Type=Integer,Description="Total allele coverage used by phasing script">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype; first allele corresponds to H1">\n')
        out.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set; all variants set to 0">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        for r in df.itertuples(index=False):
            chrom = getattr(r, "CONTIG", None) if hasattr(r, "CONTIG") else None
            # min_df currently does not include contig; default to chrX for this workflow.
            chrom = chrom or "chrX"
            pos = int(r.POS)
            ref = str(r.REF)
            alt = str(r.ALT)
            h1 = str(r.H1)
            h2 = str(r.H2)
            if h1 == ref and h2 == alt:
                gt = "0|1"
            elif h1 == alt and h2 == ref:
                gt = "1|0"
            else:
                # Should not happen for biallelic SNPs, but keep record with missing GT.
                gt = ".|."

            gene = str(r.Gene) if pd.notnull(r.Gene) else "."
            method = str(r.Stage2_method) if pd.notnull(r.Stage2_method) and str(r.Stage2_method) else "."
            stage = str(r.Stage) if pd.notnull(r.Stage) else "."
            cov = int(r.Total_coverage) if pd.notnull(r.Total_coverage) else 0
            info = f"GENE={gene};STAGE={stage};STAGE2_METHOD={method};TOTAL_COVERAGE={cov}"
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:PS\t{gt}:0\n")


# -----------------------
# Stage-2 scoring
# -----------------------
def compute_stage2_orientation_scores(
    original_filtered_formatted_df: pd.DataFrame,
    discarded_positions: List[int],
    active_calls: pd.Series,
    mode: str = "cell",
) -> pd.DataFrame:
    """
    Score both possible Stage-2 orientations.

    Orientation A:
        ALT -> H1/X1
        REF -> H2/X2

    Orientation B:
        REF -> H1/X1
        ALT -> H2/X2

    Returns:
        best_orientation:
            "alt_h1" or "ref_h1"
    """
    active_calls = active_calls.astype(str).reindex(original_filtered_formatted_df.columns)
    rows = []
    have = set(original_filtered_formatted_df.index)

    for pos in sorted(set(int(p) for p in discarded_positions)):
        ref_key = f"{pos}_ref"
        alt_key = f"{pos}_alt"

        if ref_key not in have or alt_key not in have:
            rows.append({
                "POS": pos,
                "best_orientation": "unresolved",
                "score_alt_h1": np.nan,
                "score_ref_h1": np.nan,
                "evidence": 0.0,
                "margin": np.nan,
                "n_informative_cells_stage2": 0,
                "stage2_mode": mode,
            })
            continue

        ref_counts = original_filtered_formatted_df.loc[ref_key].astype(float)
        alt_counts = original_filtered_formatted_df.loc[alt_key].astype(float)

        valid_active = active_calls.isin(["X1", "X2"])
        non_tie = ref_counts.ne(alt_counts)
        informative = valid_active & non_tie

        if informative.sum() == 0:
            rows.append({
                "POS": pos,
                "best_orientation": "unresolved",
                "score_alt_h1": np.nan,
                "score_ref_h1": np.nan,
                "evidence": 0.0,
                "margin": np.nan,
                "n_informative_cells_stage2": 0,
                "stage2_mode": mode,
            })
            continue

        act = active_calls[informative]
        ref = ref_counts[informative]
        alt = alt_counts[informative]

        if mode == "cell":
            weights = pd.Series(1.0, index=act.index)
        elif mode == "expression":
            weights = (ref - alt).abs().astype(float)
        else:
            raise ValueError(f"Unknown stage2 mode: {mode}")

        # Orientation A: ALT is H1/X1, REF is H2/X2
        # Concordant when active X1 has ALT>REF, or active X2 has REF>ALT.
        alt_h1_concordant = (
            ((act == "X1") & (alt > ref)) |
            ((act == "X2") & (ref > alt))
        )

        # Orientation B: REF is H1/X1, ALT is H2/X2
        # Concordant when active X1 has REF>ALT, or active X2 has ALT>REF.
        ref_h1_concordant = (
            ((act == "X1") & (ref > alt)) |
            ((act == "X2") & (alt > ref))
        )

        score_alt_h1 = float(weights[alt_h1_concordant].sum())
        score_ref_h1 = float(weights[ref_h1_concordant].sum())
        total = score_alt_h1 + score_ref_h1

        if total <= 0:
            best = "unresolved"
            margin = np.nan
        elif score_alt_h1 > score_ref_h1:
            best = "alt_h1"
            margin = (score_alt_h1 - score_ref_h1) / total
        elif score_ref_h1 > score_alt_h1:
            best = "ref_h1"
            margin = (score_ref_h1 - score_alt_h1) / total
        else:
            best = "tie"
            margin = 0.0

        rows.append({
            "POS": pos,
            "best_orientation": best,
            "score_alt_h1": score_alt_h1,
            "score_ref_h1": score_ref_h1,
            "evidence": total,
            "margin": margin,
            "n_informative_cells_stage2": int(informative.sum()),
            "stage2_mode": mode,
        })

    return pd.DataFrame(rows).set_index("POS")

# -----------------------
# Main
# -----------------------
def main(argv=None):
    args = parse_args(argv)
    filepath = args.filepath
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    lower_cutoff = float(args.lower_cutoff)
    upper_cutoff = 1.0 - lower_cutoff if args.upper_cutoff is None else float(args.upper_cutoff)
    if not (0.0 <= lower_cutoff <= 1.0 and 0.0 <= upper_cutoff <= 1.0):
        raise ValueError("Cutoffs must be between 0 and 1.")
    if lower_cutoff > upper_cutoff:
        raise ValueError("lower_cutoff must be <= upper_cutoff.")

    min_reads = int(args.min_reads)
    lc_tag = fmt_lc(lower_cutoff)

    if args.edge_weight_mode is None:
        edge_weight_mode = "cosine" if args.partition_mode == "expression" else "binary"
    else:
        edge_weight_mode = args.edge_weight_mode
    use_partition_weights = args.partition_mode != "unweighted"

    print(f"Using partition_mode={args.partition_mode}; edge_weight_mode={edge_weight_mode}; use_partition_weights={use_partition_weights}")

    df_std = load_and_standardize(filepath)
    debug_duplicates(df_std, outdir, examples_per_cause=10)

    pos_gene_dict = build_pos_gene_dict(df_std)
    print(f"Built support-aware POS->Gene map for {len(pos_gene_dict):,} positions")

    # ---- Stage-1 connectivity check
    run_stage1_connectivity_check(
        df_std, pos_gene_dict, min_reads, lower_cutoff, upper_cutoff, outdir, lc_tag, edge_weight_mode=edge_weight_mode
    )

    # path for Stage-1 partition/removal log
    partition_log_path = os.path.join(outdir, f"stage1_partition_log_min{min_reads}_{lc_tag}.tsv")

    # ---- Stage 1 phasing (no repair / balancing)
    (
        filtered_formatted_df,
        haplotype_1,
        haplotype_2,
        g,
        partition_array,
        all_genes,
        original_filtered_formatted_df,
    ) = generate_haplotype(
        df_std,
        pos_gene_dict,
        min_reads=min_reads,
        lower_cutoff=lower_cutoff,
        upper_cutoff=upper_cutoff,
        log_path=partition_log_path,
        edge_weight_mode=edge_weight_mode,
        use_partition_weights=use_partition_weights,
    )

    # ---- Active X calls per cell
    h1c, h2c, haplotype_sums_df, active_X_summary, c1, c2 = calculate_Xi(
        filepath, filtered_formatted_df, haplotype_1, haplotype_2, pos_gene_dict
    )
    haplotype_sums_df.to_csv(
        os.path.join(outdir, f"haplotype_sums_df_min{min_reads}_{lc_tag}.csv")
    )

    # ---- Stage 2 fill-in via concordance vs per-cell active_X
    original_snps = original_filtered_formatted_df.index.to_list()
    discarded_positions = sorted({
        int(s.split("_")[0])
        for s in original_snps
        if s not in haplotype_1 and s not in haplotype_2
    })

    stage2_scores_df = compute_stage2_orientation_scores(
        original_filtered_formatted_df=original_filtered_formatted_df,
        discarded_positions=discarded_positions,
        active_calls=haplotype_sums_df["active_X"],
        mode=args.stage2_mode,
    )

    stage2_scores_out = os.path.join(outdir, f"stage2_scores_min{min_reads}_{lc_tag}.tsv")
    stage2_scores_df.to_csv(stage2_scores_out, sep="\t")
    print(f"[DONE] Wrote {stage2_scores_out}")

    stage_map = {}
    stage1_positions = set(int(x.split("_")[0]) for x in (list(haplotype_1) + list(haplotype_2)))
    for p in stage1_positions:
        stage_map[p] = {"stage": 1, "method": None}

    haplotype_1_all = copy.deepcopy(haplotype_1)
    haplotype_2_all = copy.deepcopy(haplotype_2)

    for position, row in stage2_scores_df.iterrows():
        pos_int = int(position)
        pos_alt = f"{pos_int}_alt"
        pos_ref = f"{pos_int}_ref"

        best = row["best_orientation"]
        evidence = float(row.get("evidence", 0.0))
        margin = row.get("margin", np.nan)

        if pd.isna(margin):
            margin = 0.0

        if evidence < float(args.stage2_min_evidence):
            stage_map[pos_int] = {"stage": 2, "method": f"unresolved_low_evidence_{args.stage2_mode}"}
            continue

        if best == "alt_h1":
            haplotype_1_all.append(pos_alt)
            haplotype_2_all.append(pos_ref)
            stage_map[pos_int] = {"stage": 2, "method": f"calculated_{args.stage2_mode}_alt_h1"}

        elif best == "ref_h1":
            haplotype_1_all.append(pos_ref)
            haplotype_2_all.append(pos_alt)
            stage_map[pos_int] = {"stage": 2, "method": f"calculated_{args.stage2_mode}_ref_h1"}

        elif best == "tie":
            if args.stage2_tie_action == "skip":
                stage_map[pos_int] = {"stage": 2, "method": f"unresolved_tie_{args.stage2_mode}"}
                continue

            if random.random() < 0.5:
                haplotype_1_all.append(pos_alt)
                haplotype_2_all.append(pos_ref)
                rand_orientation = "alt_h1"
            else:
                haplotype_1_all.append(pos_ref)
                haplotype_2_all.append(pos_alt)
                rand_orientation = "ref_h1"

            stage_map[pos_int] = {"stage": 2, "method": f"random_{args.stage2_mode}_{rand_orientation}"}

        else:
            stage_map[pos_int] = {"stage": 2, "method": f"unresolved_{args.stage2_mode}"}

    # ---- Per-position total coverage
    row_totals = original_filtered_formatted_df.sum(axis=1)
    pos_series = original_filtered_formatted_df.index.to_series().str.split("_").str[0].astype(int)
    coverage_map = row_totals.groupby(pos_series).sum().to_dict()

    # ---- Global Xa/Xi labels (not saved, but computed)
    active_calls = haplotype_sums_df["active_X"]
    x1_count = int((active_calls == "X1").sum())
    x2_count = int((active_calls == "X2").sum())
    xa_h = "H1" if x1_count >= x2_count else "H2"
    xi_h = "H2" if xa_h == "H1" else "H1"
    _ = (xa_h, xi_h)

    # ---- per-SNP Xa/Xi totals
    xa_xi_by_pos = compute_xa_xi_reads(
        original_filtered_formatted_df, haplotype_1_all, haplotype_2_all, active_calls
    )

    # ---- Validation matrix
    validation_df = assemble_validation_matrix(
        filtered_formatted_df=filtered_formatted_df,
        original_filtered_formatted_df=original_filtered_formatted_df,
        g=g,
        part=partition_array,
        haplotype_1_all=haplotype_1_all,
        haplotype_2_all=haplotype_2_all,
        active_calls=active_calls,
        xa_xi_reads=xa_xi_by_pos,
    )
    validation_out = os.path.join(outdir, f"validation_metrics_min{min_reads}_{lc_tag}.tsv")
    validation_df.to_csv(validation_out, sep="\t")

    # ---- Minimal per-SNP table
    min_df = build_min_haplotype_table(
        df_std=df_std,
        haplotype_1_all=haplotype_1_all,
        haplotype_2_all=haplotype_2_all,
        pos_gene_dict=pos_gene_dict,
        stage_map=stage_map,
        coverage_map=coverage_map,
    )
    out_path = os.path.join(outdir, f"haplotypes_min{min_reads}_{lc_tag}.csv")
    min_df.to_csv(out_path, index=False)

    vcf_out = args.vcf_out or os.path.join(outdir, f"haplotypes_min{min_reads}_{lc_tag}.vcf")
    write_phased_vcf(min_df, vcf_out, sample_name=args.sample_name)
    # bgzip + tabix index VCF
    if not vcf_out.endswith(".gz"):
        subprocess.run(["bgzip", "-f", vcf_out], check=True)
        vcf_out_gz = vcf_out + ".gz"
    else:
        vcf_out_gz = vcf_out

    subprocess.run(["bcftools","index", "-t", vcf_out_gz], check=True)

    print(f"[DONE] Wrote indexed VCF: {vcf_out_gz}")
    print(f"[DONE] Wrote index: {vcf_out_gz}.tbi")
    print(f"[DONE] Wrote {out_path}")
    print(f"[DONE] Wrote {validation_out}")

    print("Done.")

if __name__ == "__main__":
    raise SystemExit(main())
