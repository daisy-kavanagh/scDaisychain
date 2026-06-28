#!/usr/bin/env python3
"""
Xi downstream pipeline in Python.

This rewrites the R-style downstream stages:

00_gene_celltype_donor_atlas.tsv
  from AnnData layers containing Xi/Xa allele-assigned counts

01_master_chrX.tsv
  by joining atlas rows to an annotation table containing PAR / escape classes

02_row_level_corrected.tsv
  by filtering rows and applying beta-prior shrinkage to Xi fractions

02_pooled_backbone_summary.tsv
  by summarising Xi recovery per donor x cell type and pooling noisy donor
  estimates toward a disease-group/cell-type backbone

03_gene_celltype_escape_summary.tsv and 03_gene_escape_pattern_summary.tsv
  by applying donor-supported or mean-based escape calls.

Expected AnnData:
  - cells x genes
  - adata.layers[Xi] and adata.layers[Xa] contain allele-assigned counts
  - adata.obs contains donor and cell-type columns
  - adata.var or adata.var_names contains gene symbols

Example:
  python xi_downstream_from_anndata.py \
    --h5ad mouse_xi_counts.h5ad \
    --annotation mouse_mouse_chrX_PAR_escape_annotation_for_xi_pipeline.tsv \
    --outdir xi_outputs \
    --xi-layer Xi \
    --xa-layer Xa \
    --donor-col donor_id \
    --celltype-col cell_type \
    --gene-name-col gene_name \
    --healthy-only false \
    --chrX-only true
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from scipy import sparse


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def str_to_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    x = str(x).strip().lower()
    if x in {"true", "t", "1", "yes", "y"}:
        return True
    if x in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got {x!r}")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False)


def as_1d_array(x) -> np.ndarray:
    """Convert sparse/dense matrix result to a flat numpy array."""
    if sparse.issparse(x):
        return np.asarray(x).ravel()
    return np.asarray(x).ravel()


def clean_chr(x: pd.Series) -> pd.Series:
    return x.astype(str).str.replace("^chr", "", regex=True)


def mean_no_na(x: pd.Series | np.ndarray) -> float:
    arr = pd.to_numeric(pd.Series(x), errors="coerce").dropna()
    if arr.empty:
        return np.nan
    return float(arr.mean())


def unique_join(x: Iterable) -> str:
    vals = sorted({str(v) for v in x if pd.notna(v)})
    return ",".join(vals)


def coalesce_series(*series: pd.Series) -> pd.Series:
    if not series:
        raise ValueError("coalesce_series requires at least one series")
    out = series[0].copy()
    for s in series[1:]:
        out = out.where(out.notna(), s)
    return out


# -----------------------------------------------------------------------------
# Beta-prior shrinkage
# -----------------------------------------------------------------------------

def estimate_beta_prior_moments(x: pd.Series | np.ndarray,
                                fallback_strength: float = 20.0,
                                min_param: float = 1e-3) -> dict[str, float]:
    """
    Estimate Beta(alpha, beta) by method of moments.

    For mean m and variance v:
      common = m(1-m)/v - 1
      alpha = m * common
      beta  = (1-m) * common

    If the empirical variance is incompatible with a beta distribution
    or is essentially zero, fall back to a beta prior centred at the
    empirical mean with total strength fallback_strength.
    """
    arr = pd.to_numeric(pd.Series(x), errors="coerce").dropna().to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    arr = arr[(arr >= 0) & (arr <= 1)]

    if arr.size == 0:
        m = 0.1
        alpha = m * fallback_strength
        beta = (1 - m) * fallback_strength
        return {
            "alpha": float(max(alpha, min_param)),
            "beta": float(max(beta, min_param)),
            "mean": float(m),
            "variance": np.nan,
            "method": "fallback_empty",
        }

    m = float(np.mean(arr))
    m = min(max(m, min_param), 1 - min_param)

    v = float(np.var(arr, ddof=1)) if arr.size > 1 else 0.0
    max_v = m * (1 - m)

    if not np.isfinite(v) or v <= 0 or v >= max_v:
        alpha = m * fallback_strength
        beta = (1 - m) * fallback_strength
        method = "fallback_strength"
    else:
        common = max_v / v - 1
        alpha = m * common
        beta = (1 - m) * common
        method = "method_of_moments"

    alpha = float(max(alpha, min_param))
    beta = float(max(beta, min_param))

    return {
        "alpha": alpha,
        "beta": beta,
        "mean": float(alpha / (alpha + beta)),
        "variance": float(v),
        "method": method,
    }


# -----------------------------------------------------------------------------
# Escape-call helpers, mirroring 00_xi_config_helpers.R
# -----------------------------------------------------------------------------

def resolve_escape_params(mode: Optional[str],
                          xi_cutoff: Optional[float],
                          min_donors: Optional[int]) -> dict:
    mode = mode or "recovery"
    if mode not in {"recovery", "strict", "mean_strict"}:
        raise ValueError(f"Unknown escape call mode: {mode}. Use recovery, strict, or mean_strict.")

    if xi_cutoff is None:
        xi_cutoff = {
            "recovery": 0.10,
            "strict": 0.20,
            "mean_strict": 0.20,
        }[mode]

    if min_donors is None:
        min_donors = {
            "recovery": 1,
            "strict": 2,
            "mean_strict": 2,
        }[mode]

    return {
        "mode": mode,
        "xi_cutoff": float(xi_cutoff),
        "min_donors": int(min_donors),
    }


def apply_escape_call_to_summary(d: pd.DataFrame,
                                 mode: str,
                                 xi_cutoff: float,
                                 min_donors: int) -> pd.DataFrame:
    """
    Apply a consistent escape call to a gene x cell-type summary table.

    Required columns:
      - mean_xi
      - n_donors

    Preferred for recovery/strict:
      - n_escape_donors
    """
    required = {"mean_xi", "n_donors"}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"Summary table needs columns: {', '.join(sorted(missing))}")

    out = d.copy()
    out["informative"] = out["n_donors"] >= min_donors

    if mode in {"recovery", "strict"}:
        if "n_escape_donors" in out.columns:
            out["escape_call"] = out["n_escape_donors"] >= min_donors
        else:
            print(
                f"WARNING: n_escape_donors absent; falling back to mean_strict-style call for mode={mode}",
                file=sys.stderr,
            )
            out["escape_call"] = out["informative"] & (out["mean_xi"] >= xi_cutoff)
    elif mode == "mean_strict":
        out["escape_call"] = out["informative"] & (out["mean_xi"] >= xi_cutoff)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return out


def classify_escape_pattern_from_counts(n_escape_celltypes: int,
                                        n_informative_celltypes: int,
                                        n_non_escape_informative_celltypes: int,
                                        shared_min_escape_celltypes: int = 6,
                                        intermediate_min_escape_celltypes: int = 3,
                                        intermediate_max_escape_celltypes: int = 5,
                                        intermediate_min_informative_celltypes: int = 0,
                                        intermediate_min_non_escape_celltypes: int = 0,
                                        selective_max_escape_celltypes: int = 2,
                                        selective_min_informative_celltypes: int = 0,
                                        selective_min_non_escape_celltypes: int = 0,
                                        require_non_escape_for_modulated: bool = False) -> str:
    shared = n_escape_celltypes >= shared_min_escape_celltypes

    if require_non_escape_for_modulated:
        intermediate = (
            n_escape_celltypes >= intermediate_min_escape_celltypes
            and n_escape_celltypes <= intermediate_max_escape_celltypes
            and n_informative_celltypes >= intermediate_min_informative_celltypes
            and n_non_escape_informative_celltypes >= intermediate_min_non_escape_celltypes
        )
        selective = (
            n_escape_celltypes >= 1
            and n_escape_celltypes <= selective_max_escape_celltypes
            and n_informative_celltypes >= selective_min_informative_celltypes
            and n_non_escape_informative_celltypes >= selective_min_non_escape_celltypes
        )
    else:
        intermediate = (
            n_escape_celltypes >= intermediate_min_escape_celltypes
            and n_escape_celltypes < shared_min_escape_celltypes
        )
        selective = (
            n_escape_celltypes >= 1
            and n_escape_celltypes < intermediate_min_escape_celltypes
        )

    if shared:
        return "shared_constitutive"
    if intermediate:
        return "intermediate_strict"
    if selective:
        return "selective_strict"
    return "does_not_pass_strict"


# -----------------------------------------------------------------------------
# Stage 00: AnnData -> gene x cell-type x donor atlas
# -----------------------------------------------------------------------------

def load_anndata(h5ad: Path):
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("Install anndata first: pip install anndata") from exc
    return ad.read_h5ad(h5ad)


def get_layer(adata, layer_name: str):
    if layer_name not in adata.layers:
        raise KeyError(
            f"AnnData layer {layer_name!r} not found. Available layers: {list(adata.layers.keys())}"
        )
    return adata.layers[layer_name]


def get_gene_names(adata, gene_name_col: Optional[str]) -> pd.Series:
    """
    Return gene symbols/names aligned to adata.var rows.

    Use adata.var_names when:
      - --gene-name-col is omitted
      - --gene-name-col is one of: index, var_names, var_index
      - --gene-name-col was set to gene_name but adata.var has no gene_name column

    The last fallback is intentional because Cell Ranger-derived AnnData objects
    often store gene symbols in var_names and Ensembl IDs in adata.var["gene_ids"].
    """
    use_index_tokens = {"index", "var_names", "var_index", "adata.var_names"}

    if gene_name_col is None or str(gene_name_col) in use_index_tokens:
        genes = pd.Series(adata.var_names.astype(str), index=adata.var_names)
    elif gene_name_col in adata.var.columns:
        genes = adata.var[gene_name_col].astype(str)
    else:
        available = list(adata.var.columns)

        if gene_name_col == "gene_name":
            print(
                "WARNING: --gene-name-col gene_name was requested, but adata.var has no "
                "gene_name column. Falling back to adata.var_names. Available adata.var "
                f"columns: {available}",
                file=sys.stderr,
            )
            genes = pd.Series(adata.var_names.astype(str), index=adata.var_names)
        else:
            raise KeyError(
                f"gene_name_col {gene_name_col!r} not found in adata.var. "
                f"Available columns: {available}. To use the AnnData index, omit "
                "--gene-name-col or set --gene-name-col var_names."
            )

    if genes.isna().any():
        raise ValueError("Some gene names are NA.")
    if genes.duplicated().any():
        duplicated = genes[genes.duplicated()].head(10).tolist()
        print(
            "WARNING: gene names are duplicated. Output is per var row; downstream merge by gene_name "
            f"may collapse/duplicate these. Examples: {duplicated}",
            file=sys.stderr,
        )
    return genes.reset_index(drop=True)


def build_gene_celltype_donor_atlas(
    adata,
    xi_layer: str,
    xa_layer: str,
    donor_col: str,
    celltype_col: str,
    gene_name_col: Optional[str],
    disease_group_col: Optional[str],
    sample_group_col: Optional[str],
    age_group_col: Optional[str],
    min_total_reads: int,
    min_cells_with_signal: int,
    keep_zero_rows: bool = False,
) -> pd.DataFrame:
    Xi = get_layer(adata, xi_layer)
    Xa = get_layer(adata, xa_layer)

    if Xi.shape != adata.shape or Xa.shape != adata.shape:
        raise ValueError(
            f"Layer shapes must match adata.shape={adata.shape}; "
            f"got Xi={Xi.shape}, Xa={Xa.shape}"
        )

    if not sparse.issparse(Xi):
        Xi = sparse.csr_matrix(Xi)
    else:
        Xi = Xi.tocsr()

    if not sparse.issparse(Xa):
        Xa = sparse.csr_matrix(Xa)
    else:
        Xa = Xa.tocsr()

    obs = adata.obs.copy()
    required_obs = [donor_col, celltype_col]
    for col in required_obs:
        if col not in obs.columns:
            raise KeyError(f"Required adata.obs column not found: {col}")

    # Optional columns are retained if present; otherwise filled.
    optional_map = {
        "disease_group": disease_group_col,
        "sample_group": sample_group_col,
        "age_group": age_group_col,
    }
    for canonical, col in optional_map.items():
        if col and col in obs.columns:
            obs[canonical] = obs[col].astype("string")
        else:
            obs[canonical] = pd.NA

    obs["_donor_id"] = obs[donor_col].astype("string")
    obs["_cell_type"] = obs[celltype_col].astype("string")

    genes = get_gene_names(adata, gene_name_col)
    gene_df = pd.DataFrame({
        "var_index": np.arange(adata.n_vars, dtype=int),
        "gene_name": genes.to_numpy(),
    })

    group_cols = ["_donor_id", "_cell_type", "sample_group", "disease_group", "age_group"]
    # dropna=False is important: keep groups even if optional metadata are NA.
    groups = obs.reset_index(drop=True).groupby(group_cols, dropna=False).indices

    rows = []
    for group_key, cell_idx in groups.items():
        donor_id, cell_type, sample_group, disease_group, age_group = group_key
        cell_idx = np.asarray(cell_idx, dtype=int)

        Xi_g = Xi[cell_idx, :]
        Xa_g = Xa[cell_idx, :]
        total_g = Xi_g + Xa_g

        xi_reads = as_1d_array(Xi_g.sum(axis=0))
        xa_reads = as_1d_array(Xa_g.sum(axis=0))
        total_reads = xi_reads + xa_reads
        n_cells_with_signal = np.asarray(total_g.getnnz(axis=0)).ravel().astype(int)

        if keep_zero_rows:
            keep = np.ones(adata.n_vars, dtype=bool)
        else:
            keep = total_reads > 0

        if not np.any(keep):
            continue

        out = gene_df.loc[keep, ["gene_name"]].copy()
        out["donor_id"] = str(donor_id) if pd.notna(donor_id) else pd.NA
        out["cell_type"] = str(cell_type) if pd.notna(cell_type) else pd.NA
        out["sample_group"] = str(sample_group) if pd.notna(sample_group) else pd.NA
        out["disease_group"] = str(disease_group) if pd.notna(disease_group) else pd.NA
        out["age_group"] = str(age_group) if pd.notna(age_group) else pd.NA
        out["n_cells_total_group"] = int(len(cell_idx))
        out["xi_reads"] = xi_reads[keep].astype(float)
        out["xa_reads"] = xa_reads[keep].astype(float)
        out["total_allelic_reads"] = total_reads[keep].astype(float)
        out["n_cells_with_signal"] = n_cells_with_signal[keep].astype(int)
        out["informative_firstpass"] = (
            (out["total_allelic_reads"] >= min_total_reads)
            & (out["n_cells_with_signal"] >= min_cells_with_signal)
        )
        rows.append(out)

    if not rows:
        return pd.DataFrame(
            columns=[
                "gene_name", "donor_id", "cell_type", "sample_group", "disease_group",
                "age_group", "n_cells_total_group", "xi_reads", "xa_reads",
                "total_allelic_reads", "n_cells_with_signal", "informative_firstpass",
            ]
        )

    atlas = pd.concat(rows, ignore_index=True)
    atlas["xi_raw_firstpass"] = np.where(
        atlas["total_allelic_reads"] > 0,
        atlas["xi_reads"] / atlas["total_allelic_reads"],
        np.nan,
    )
    return atlas


# -----------------------------------------------------------------------------
# Stage 01: annotate and build master chrX table
# -----------------------------------------------------------------------------

def map_escape_display(x: pd.Series) -> pd.Series:
    display_map = {
        "nonPAR_escape": "Escape",
        "nonPAR_variable_escape": "Variable",
        "nonPAR_subject": "Inactive",
        "nonPAR_uncertain": "Unknown",
        "PAR": "PAR",
    }
    return x.map(display_map).fillna(x)


def standardise_annotation(ann: pd.DataFrame) -> pd.DataFrame:
    ann = ann.copy()

    if "gene_symbol" in ann.columns and "gene_name" not in ann.columns:
        ann = ann.rename(columns={"gene_symbol": "gene_name"})
    if "gene_name" not in ann.columns:
        raise ValueError("Annotation table must contain gene_name or gene_symbol.")

    if "chrom" in ann.columns and "chromosome_clean" not in ann.columns:
        ann["chromosome_clean"] = clean_chr(ann["chrom"])

    if "is_chrX" not in ann.columns:
        if "chromosome_clean" in ann.columns:
            ann["is_chrX"] = ann["chromosome_clean"].astype(str).isin(["X", "x"])
        elif "chrom" in ann.columns:
            ann["is_chrX"] = ann["chrom"].astype(str).isin(["chrX", "X", "x"])

    if "preferred_escape_class" not in ann.columns:
        if "xci_status" in ann.columns:
            par_flag = pd.Series(False, index=ann.index)
            if "PAR_status" in ann.columns:
                par_flag |= ann["PAR_status"].isin(["PAR", "PAR1", "PAR2"])
            if "PAR_region" in ann.columns:
                par_flag |= ann["PAR_region"].isin(["PAR1", "PAR2"])
            par_flag |= ann["xci_status"].eq("PAR")

            ann["preferred_escape_class"] = pd.NA
            ann.loc[par_flag, "preferred_escape_class"] = "PAR"
            ann.loc[ann["xci_status"].eq("Escape"), "preferred_escape_class"] = "nonPAR_escape"
            ann.loc[ann["xci_status"].eq("Variable"), "preferred_escape_class"] = "nonPAR_variable_escape"
            ann.loc[ann["xci_status"].eq("Inactive"), "preferred_escape_class"] = "nonPAR_subject"
            ann.loc[ann["xci_status"].eq("Unknown"), "preferred_escape_class"] = "nonPAR_uncertain"

        elif "category" in ann.columns:
            par_flag = pd.Series(False, index=ann.index)
            if "PAR_status" in ann.columns:
                par_flag |= ann["PAR_status"].isin(["PAR", "PAR1", "PAR2"])
            if "PAR_region" in ann.columns:
                par_flag |= ann["PAR_region"].isin(["PAR1", "PAR2"])
            par_flag |= ann["category"].eq("PAR")

            ann["preferred_escape_class"] = pd.NA
            ann.loc[par_flag, "preferred_escape_class"] = "PAR"
            ann.loc[ann["category"].eq("Escape"), "preferred_escape_class"] = "nonPAR_escape"
            ann.loc[ann["category"].eq("Variable"), "preferred_escape_class"] = "nonPAR_variable_escape"
            ann.loc[ann["category"].eq("Inactive"), "preferred_escape_class"] = "nonPAR_subject"
            ann.loc[ann["category"].eq("Unknown"), "preferred_escape_class"] = "nonPAR_uncertain"

    if "atlas_gene_group" not in ann.columns and "preferred_escape_class" in ann.columns:
        ann["atlas_gene_group"] = ann["preferred_escape_class"]

    # Avoid duplicate annotation rows silently expanding the master table.
    # Keep the first but write duplicate info in QC outside this function.
    ann = ann.drop_duplicates(subset=["gene_name"], keep="first")
    return ann


def annotate_and_build_master(
    atlas: pd.DataFrame,
    annotation_path: Path,
    processed_dir: Path,
    tables_dir: Path,
    healthy_only: bool,
    chrX_only: bool,
) -> pd.DataFrame:
    ann_raw = pd.read_csv(annotation_path, sep=None, engine="python", dtype=str)
    duplicate_genes = ann_raw["gene_name"].duplicated().sum() if "gene_name" in ann_raw.columns else np.nan
    ann = standardise_annotation(ann_raw)

    dt = atlas.merge(ann, on="gene_name", how="left", suffixes=("", ".ann"))

    if "atlas_gene_group" in dt.columns:
        dt["atlas_gene_group_display"] = map_escape_display(dt["atlas_gene_group"])

    save_tsv(dt, processed_dir / "01_master_pre_filter.tsv")

    if "disease_group" not in dt.columns:
        print(
            "WARNING: disease_group column missing from atlas; setting to NA. "
            "healthy_only filtering will not work until atlas supplies disease_group.",
            file=sys.stderr,
        )
        dt["disease_group"] = pd.NA

    if healthy_only:
        dt = dt.loc[dt["disease_group"].eq("Healthy")].copy()

    if chrX_only:
        if "is_chrX" in dt.columns:
            # Handles bools, 1/0, and strings.
            is_chrX = (
                dt["is_chrX"].eq(True)
                | dt["is_chrX"].eq(1)
                | dt["is_chrX"].astype(str).str.lower().isin(["true", "1", "x", "chrx"])
            )
            dt = dt.loc[is_chrX].copy()
        elif "chromosome_clean" in dt.columns:
            dt = dt.loc[dt["chromosome_clean"].astype(str).isin(["X", "x", "chrX"])].copy()

    save_tsv(dt, processed_dir / "01_master_chrX.tsv")
    n_missing_pref = int(dt["preferred_escape_class"].isna().sum()) if "preferred_escape_class" in dt.columns else np.nan
    merge_qc = pd.DataFrame([{
        "annotation_path": str(annotation_path),
        "healthy_only_filter": healthy_only,
        "chrX_only_filter": chrX_only,
        "n_rows": len(dt),
        "n_unique_genes": dt["gene_name"].nunique(dropna=True),
        "n_unique_donors": dt["donor_id"].nunique(dropna=True) if "donor_id" in dt.columns else np.nan,
        "donors": unique_join(dt["donor_id"]) if "donor_id" in dt.columns else "",
        "disease_groups": unique_join(dt["disease_group"]) if "disease_group" in dt.columns else "",
        "n_duplicate_annotation_gene_names_input": duplicate_genes,
        "n_missing_atlas_gene_group": int(dt["atlas_gene_group"].isna().sum()) if "atlas_gene_group" in dt.columns else np.nan,
        "n_missing_preferred_escape_class": n_missing_pref,
    }])
    save_tsv(merge_qc, tables_dir / "01_annotation_merge_qc.tsv")

    if "sample_group" in dt.columns:
        grp_qc = (
            dt.groupby(["sample_group", "disease_group"], dropna=False)
            .agg(
                n_rows=("gene_name", "size"),
                n_donors=("donor_id", lambda x: x.nunique(dropna=True)),
                donors=("donor_id", unique_join),
            )
            .reset_index()
        )
        save_tsv(grp_qc, tables_dir / "01_sample_group_qc.tsv")

    return dt


# -----------------------------------------------------------------------------
# Stage 02: row filtering, shrinkage, pooling
# -----------------------------------------------------------------------------

def shrinkage_pooling_master_chrX(
    master: pd.DataFrame,
    processed_dir: Path,
    tables_dir: Path,
    min_total_reads: int,
    min_cells_with_signal: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    req_cols = {
        "donor_id", "cell_type", "gene_name", "xi_reads", "xa_reads",
        "total_allelic_reads", "informative_firstpass", "n_cells_with_signal",
    }
    missing = req_cols - set(master.columns)
    if missing:
        raise ValueError(f"Master table missing columns: {', '.join(sorted(missing))}")

    dt = master.copy()
    for col in ["disease_group", "sample_group", "age_group"]:
        if col not in dt.columns:
            dt[col] = pd.NA

    dt["informative_firstpass"] = dt["informative_firstpass"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
    for col in ["xi_reads", "xa_reads", "total_allelic_reads", "n_cells_with_signal"]:
        dt[col] = pd.to_numeric(dt[col], errors="coerce")

    dt = dt.loc[
        dt["informative_firstpass"]
        & (dt["total_allelic_reads"] >= min_total_reads)
        & (dt["n_cells_with_signal"] >= min_cells_with_signal)
    ].copy()

    dt["passed_row_filter"] = True
    dt["xi_raw"] = dt["xi_reads"] / dt["total_allelic_reads"]

    prior = estimate_beta_prior_moments(dt["xi_raw"])
    dt["xi_shrunk"] = (
        (dt["xi_reads"] + prior["alpha"])
        / (dt["total_allelic_reads"] + prior["alpha"] + prior["beta"])
    )

    save_tsv(pd.DataFrame([prior]), tables_dir / "02_beta_prior.tsv")
    save_tsv(dt, processed_dir / "02_row_level_corrected.tsv")

    filter_audit = pd.DataFrame([{
        "min_total_reads": min_total_reads,
        "min_cells_with_signal": min_cells_with_signal,
        "n_rows_retained": len(dt),
        "n_genes_retained": dt["gene_name"].nunique(dropna=True),
        "n_donors_retained": dt["donor_id"].nunique(dropna=True),
        "n_celltypes_retained": dt["cell_type"].nunique(dropna=True),
        "disease_groups": unique_join(dt["disease_group"]),
    }])
    save_tsv(filter_audit, tables_dir / "02_row_filter_audit.tsv")

    # Disease-specific backbone prevents one disease group from influencing the other.
    celltype_backbone = (
        dt.groupby(["disease_group", "cell_type"], dropna=False)
        .agg(
            xi_shrunk_celltype_mean=("xi_shrunk", mean_no_na),
            celltype_total_reads=("total_allelic_reads", "sum"),
        )
        .reset_index()
    )

    donor_cell = (
        dt.groupby(["donor_id", "sample_group", "disease_group", "age_group", "cell_type"], dropna=False)
        .agg(
            xi_shrunk_mean=("xi_shrunk", mean_no_na),
            n_informative_genes=("gene_name", lambda x: int(x.nunique(dropna=True))),
            total_reads=("total_allelic_reads", "sum"),
        )
        .reset_index()
    )

    donor_cell = donor_cell.merge(
        celltype_backbone,
        on=["disease_group", "cell_type"],
        how="left",
    )

    median_reads = float(np.nanmedian(donor_cell["total_reads"].to_numpy(dtype=float))) if len(donor_cell) else np.nan
    if not np.isfinite(median_reads) or median_reads <= 0:
        median_reads = 1.0

    donor_cell["pool_weight"] = donor_cell["total_reads"] / (donor_cell["total_reads"] + median_reads)
    donor_cell["xi_pooled"] = (
        donor_cell["pool_weight"] * donor_cell["xi_shrunk_mean"]
        + (1 - donor_cell["pool_weight"]) * donor_cell["xi_shrunk_celltype_mean"]
    )

    save_tsv(donor_cell, tables_dir / "02_pooled_backbone_summary.tsv")

    return dt, donor_cell, prior


# -----------------------------------------------------------------------------
# Stage 03: gene x cell-type escape summary and pattern classification
# -----------------------------------------------------------------------------

def build_escape_summaries(
    row_level: pd.DataFrame,
    processed_dir: Path,
    tables_dir: Path,
    mode: str,
    xi_cutoff: float,
    min_donors: int,
    group_by_disease: bool = True,
    shared_min_escape_celltypes: int = 6,
    intermediate_min_escape_celltypes: int = 3,
    intermediate_max_escape_celltypes: int = 5,
    selective_max_escape_celltypes: int = 2,
    require_non_escape_for_modulated: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dt = row_level.copy()

    dt["donor_escape"] = dt["xi_shrunk"] >= xi_cutoff

    annot_cols = [
        c for c in [
            "preferred_escape_class",
            "atlas_gene_group",
            "atlas_gene_group_display",
            "PAR_status",
            "chromosome_clean",
            "is_chrX",
        ]
        if c in dt.columns
    ]

    by_cols = ["gene_name", "cell_type"]
    if group_by_disease and "disease_group" in dt.columns:
        by_cols = ["disease_group"] + by_cols

    # Preserve one annotation value per gene if present.
    agg_dict = {
        "mean_xi": ("xi_shrunk", mean_no_na),
        "mean_xi_raw": ("xi_raw", mean_no_na),
        "n_donors": ("donor_id", lambda x: int(x.nunique(dropna=True))),
        "n_escape_donors": ("donor_escape", lambda x: int(np.nansum(x.astype(bool)))),
        "total_allelic_reads": ("total_allelic_reads", "sum"),
        "n_rows": ("gene_name", "size"),
    }
    for col in annot_cols:
        agg_dict[col] = (col, lambda x: x.dropna().iloc[0] if x.dropna().size else pd.NA)

    summary = dt.groupby(by_cols, dropna=False).agg(**agg_dict).reset_index()
    summary = apply_escape_call_to_summary(
        summary,
        mode=mode,
        xi_cutoff=xi_cutoff,
        min_donors=min_donors,
    )

    save_tsv(summary, processed_dir / "03_gene_celltype_escape_summary.tsv")

    pattern_by = ["gene_name"]
    if group_by_disease and "disease_group" in summary.columns:
        pattern_by = ["disease_group", "gene_name"]

    pattern_rows = []
    for keys, g in summary.groupby(pattern_by, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)

        informative = g["informative"].astype(bool)
        escape = g["escape_call"].astype(bool)

        n_escape_celltypes = int((informative & escape).sum())
        n_informative_celltypes = int(informative.sum())
        n_non_escape_informative_celltypes = int((informative & ~escape).sum())

        pattern = classify_escape_pattern_from_counts(
            n_escape_celltypes=n_escape_celltypes,
            n_informative_celltypes=n_informative_celltypes,
            n_non_escape_informative_celltypes=n_non_escape_informative_celltypes,
            shared_min_escape_celltypes=shared_min_escape_celltypes,
            intermediate_min_escape_celltypes=intermediate_min_escape_celltypes,
            intermediate_max_escape_celltypes=intermediate_max_escape_celltypes,
            selective_max_escape_celltypes=selective_max_escape_celltypes,
            require_non_escape_for_modulated=require_non_escape_for_modulated,
        )

        row = dict(zip(pattern_by, keys))
        row.update({
            "n_escape_celltypes": n_escape_celltypes,
            "n_informative_celltypes": n_informative_celltypes,
            "n_non_escape_informative_celltypes": n_non_escape_informative_celltypes,
            "escape_pattern": pattern,
            "escaping_celltypes": ",".join(sorted(g.loc[informative & escape, "cell_type"].astype(str).unique())),
            "informative_celltypes": ",".join(sorted(g.loc[informative, "cell_type"].astype(str).unique())),
        })

        # Add gene-level annotation fields if present.
        for col in annot_cols:
            vals = g[col].dropna()
            row[col] = vals.iloc[0] if len(vals) else pd.NA

        pattern_rows.append(row)

    pattern_summary = pd.DataFrame(pattern_rows)
    save_tsv(pattern_summary, processed_dir / "03_gene_escape_pattern_summary.tsv")

    call_qc = pd.DataFrame([{
        "mode": mode,
        "xi_cutoff": xi_cutoff,
        "min_donors": min_donors,
        "group_by_disease": group_by_disease,
        "n_gene_celltype_rows": len(summary),
        "n_gene_celltype_escape_calls": int(summary["escape_call"].sum()),
        "n_gene_pattern_rows": len(pattern_summary),
    }])
    save_tsv(call_qc, tables_dir / "03_escape_call_qc.tsv")

    return summary, pattern_summary


# -----------------------------------------------------------------------------
# Run parameters
# -----------------------------------------------------------------------------

def write_run_parameters(path: Path, args: argparse.Namespace, extra: Optional[dict] = None) -> None:
    vals = vars(args).copy()
    if extra:
        vals.update(extra)

    # Store compact JSON for complex objects if needed.
    rows = []
    for k, v in vals.items():
        if isinstance(v, (dict, list, tuple)):
            v = json.dumps(v)
        rows.append({"parameter": k, "value": str(v)})
    save_tsv(pd.DataFrame(rows), path)


# -----------------------------------------------------------------------------
# Callable pipeline API
# -----------------------------------------------------------------------------

def run_xi_downstream_from_adata(
    adata,
    annotation: Path | str,
    outdir: Path | str,
    xi_layer: str = "Xi",
    xa_layer: str = "Xa",
    donor_col: str = "donor_id",
    celltype_col: str = "cell_type",
    gene_name_col: Optional[str] = None,
    disease_group_col: Optional[str] = "disease_group",
    sample_group_col: Optional[str] = "sample_group",
    age_group_col: Optional[str] = "age_group",
    min_total_reads: int = 3,
    min_cells_with_signal: int = 2,
    keep_zero_rows: bool = False,
    healthy_only: bool = False,
    chrX_only: bool = True,
    escape_mode: str = "recovery",
    xi_cutoff: Optional[float] = None,
    min_donors: Optional[int] = None,
    group_escape_by_disease: bool = True,
    shared_min_escape_celltypes: int = 6,
    intermediate_min_escape_celltypes: int = 3,
    intermediate_max_escape_celltypes: int = 5,
    selective_max_escape_celltypes: int = 2,
    require_non_escape_for_modulated: bool = False,
    skip_03: bool = False,
) -> dict[str, pd.DataFrame | dict | Path]:
    """
    Run the Xi downstream pipeline from an in-memory AnnData object.

    This is the function to use inside a Jupyter notebook.

    Returns a dictionary containing the main DataFrames and paths:
      - atlas
      - master
      - row_level
      - donor_cell
      - beta_prior
      - gene_celltype_escape_summary, if skip_03 is False
      - gene_escape_pattern_summary, if skip_03 is False
      - outdir, processed_dir, tables_dir
    """
    annotation = Path(annotation)
    outdir = Path(outdir)

    processed_dir = ensure_dir(outdir / "processed")
    tables_dir = ensure_dir(outdir / "tables")

    escape_params = resolve_escape_params(
        mode=escape_mode,
        xi_cutoff=xi_cutoff,
        min_donors=min_donors,
    )

    print("Building 00_gene_celltype_donor_atlas.tsv")
    atlas = build_gene_celltype_donor_atlas(
        adata=adata,
        xi_layer=xi_layer,
        xa_layer=xa_layer,
        donor_col=donor_col,
        celltype_col=celltype_col,
        gene_name_col=gene_name_col,
        disease_group_col=disease_group_col,
        sample_group_col=sample_group_col,
        age_group_col=age_group_col,
        min_total_reads=min_total_reads,
        min_cells_with_signal=min_cells_with_signal,
        keep_zero_rows=keep_zero_rows,
    )
    save_tsv(atlas, processed_dir / "00_gene_celltype_donor_atlas.tsv")

    print("Building 01_master_chrX.tsv")
    master = annotate_and_build_master(
        atlas=atlas,
        annotation_path=annotation,
        processed_dir=processed_dir,
        tables_dir=tables_dir,
        healthy_only=healthy_only,
        chrX_only=chrX_only,
    )

    print("Building 02_row_level_corrected.tsv and 02_pooled_backbone_summary.tsv")
    row_level, donor_cell, prior = shrinkage_pooling_master_chrX(
        master=master,
        processed_dir=processed_dir,
        tables_dir=tables_dir,
        min_total_reads=min_total_reads,
        min_cells_with_signal=min_cells_with_signal,
    )

    result = {
        "atlas": atlas,
        "master": master,
        "row_level": row_level,
        "donor_cell": donor_cell,
        "beta_prior": prior,
        "outdir": outdir,
        "processed_dir": processed_dir,
        "tables_dir": tables_dir,
    }

    if not skip_03:
        print("Building 03_gene_celltype_escape_summary.tsv and 03_gene_escape_pattern_summary.tsv")
        gene_celltype_summary, gene_pattern_summary = build_escape_summaries(
            row_level=row_level,
            processed_dir=processed_dir,
            tables_dir=tables_dir,
            mode=escape_params["mode"],
            xi_cutoff=escape_params["xi_cutoff"],
            min_donors=escape_params["min_donors"],
            group_by_disease=group_escape_by_disease,
            shared_min_escape_celltypes=shared_min_escape_celltypes,
            intermediate_min_escape_celltypes=intermediate_min_escape_celltypes,
            intermediate_max_escape_celltypes=intermediate_max_escape_celltypes,
            selective_max_escape_celltypes=selective_max_escape_celltypes,
            require_non_escape_for_modulated=require_non_escape_for_modulated,
        )
        result["gene_celltype_escape_summary"] = gene_celltype_summary
        result["gene_escape_pattern_summary"] = gene_pattern_summary

    run_params = {
        "annotation": str(annotation),
        "outdir": str(outdir),
        "xi_layer": xi_layer,
        "xa_layer": xa_layer,
        "donor_col": donor_col,
        "celltype_col": celltype_col,
        "gene_name_col": gene_name_col,
        "disease_group_col": disease_group_col,
        "sample_group_col": sample_group_col,
        "age_group_col": age_group_col,
        "min_total_reads": min_total_reads,
        "min_cells_with_signal": min_cells_with_signal,
        "keep_zero_rows": keep_zero_rows,
        "healthy_only": healthy_only,
        "chrX_only": chrX_only,
        "escape_mode": escape_mode,
        "resolved_escape_mode": escape_params["mode"],
        "resolved_xi_cutoff": escape_params["xi_cutoff"],
        "resolved_min_donors": escape_params["min_donors"],
        "group_escape_by_disease": group_escape_by_disease,
        "shared_min_escape_celltypes": shared_min_escape_celltypes,
        "intermediate_min_escape_celltypes": intermediate_min_escape_celltypes,
        "intermediate_max_escape_celltypes": intermediate_max_escape_celltypes,
        "selective_max_escape_celltypes": selective_max_escape_celltypes,
        "require_non_escape_for_modulated": require_non_escape_for_modulated,
        "skip_03": skip_03,
        "beta_prior_alpha": prior["alpha"],
        "beta_prior_beta": prior["beta"],
        "beta_prior_method": prior["method"],
    }
    save_tsv(
        pd.DataFrame(
            [{"parameter": k, "value": str(v)} for k, v in run_params.items()]
        ),
        tables_dir / "run_parameters.tsv",
    )

    print()
    print("Done.")
    print(f"Output directory: {outdir}")
    print(f"Rows in atlas: {len(atlas):,}")
    print(f"Rows in master: {len(master):,}")
    print(f"Rows after 02 filtering: {len(row_level):,}")
    print(f"Donor/cell-type pooled rows: {len(donor_cell):,}")

    return result


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Xi downstream pipeline from AnnData Xi/Xa layers. "
            "Can be run as a CLI from .h5ad, or imported and called on an in-memory AnnData object."
        )
    )

    parser.add_argument("--h5ad", required=True, type=Path, help="Input AnnData .h5ad file.")
    parser.add_argument("--annotation", required=True, type=Path, help="Annotation TSV/CSV with gene_name and escape/PAR classes.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")

    parser.add_argument("--xi-layer", default="Xi", help="AnnData layer containing Xi counts. Default: Xi")
    parser.add_argument("--xa-layer", default="Xa", help="AnnData layer containing Xa counts. Default: Xa")

    parser.add_argument("--donor-col", default="donor_id", help="adata.obs donor column. Default: donor_id")
    parser.add_argument("--celltype-col", default="cell_type", help="adata.obs cell-type column. Default: cell_type")
    parser.add_argument("--gene-name-col", default=None, help="adata.var gene symbol/name column. Default: use adata.var_names. You can also set var_names/index explicitly.")

    parser.add_argument("--disease-group-col", default="disease_group", help="Optional adata.obs disease group column.")
    parser.add_argument("--sample-group-col", default="sample_group", help="Optional adata.obs sample group column.")
    parser.add_argument("--age-group-col", default="age_group", help="Optional adata.obs age group column.")

    parser.add_argument("--min-total-reads", type=int, default=3, help="Minimum total allelic reads per donor/celltype/gene row.")
    parser.add_argument("--min-cells-with-signal", type=int, default=2, help="Minimum cells with Xi+Xa signal per donor/celltype/gene row.")
    parser.add_argument("--keep-zero-rows", type=str_to_bool, default=False, help="Write zero-count gene rows in atlas. Default: false")

    parser.add_argument("--healthy-only", type=str_to_bool, default=False, help="Keep only disease_group == Healthy after annotation. Default: false")
    parser.add_argument("--chrX-only", type=str_to_bool, default=True, help="Keep only chrX genes after annotation. Default: true")

    parser.add_argument("--escape-mode", default="recovery", choices=["recovery", "strict", "mean_strict"], help="Escape call mode.")
    parser.add_argument("--xi-cutoff", type=float, default=None, help="Xi cutoff for donor escape calls. Defaults depend on mode.")
    parser.add_argument("--min-donors", type=int, default=None, help="Minimum donors for escape call. Defaults depend on mode.")
    parser.add_argument("--group-escape-by-disease", type=str_to_bool, default=True, help="Build escape summaries separately by disease_group. Default: true")

    parser.add_argument("--shared-min-escape-celltypes", type=int, default=6)
    parser.add_argument("--intermediate-min-escape-celltypes", type=int, default=3)
    parser.add_argument("--intermediate-max-escape-celltypes", type=int, default=5)
    parser.add_argument("--selective-max-escape-celltypes", type=int, default=2)
    parser.add_argument("--require-non-escape-for-modulated", type=str_to_bool, default=False)

    parser.add_argument("--skip-03", action="store_true", help="Stop after stage 02.")

    args = parser.parse_args()

    print(f"Reading AnnData: {args.h5ad}")
    adata = load_anndata(args.h5ad)

    run_xi_downstream_from_adata(
        adata=adata,
        annotation=args.annotation,
        outdir=args.outdir,
        xi_layer=args.xi_layer,
        xa_layer=args.xa_layer,
        donor_col=args.donor_col,
        celltype_col=args.celltype_col,
        gene_name_col=args.gene_name_col,
        disease_group_col=args.disease_group_col,
        sample_group_col=args.sample_group_col,
        age_group_col=args.age_group_col,
        min_total_reads=args.min_total_reads,
        min_cells_with_signal=args.min_cells_with_signal,
        keep_zero_rows=args.keep_zero_rows,
        healthy_only=args.healthy_only,
        chrX_only=args.chrX_only,
        escape_mode=args.escape_mode,
        xi_cutoff=args.xi_cutoff,
        min_donors=args.min_donors,
        group_escape_by_disease=args.group_escape_by_disease,
        shared_min_escape_celltypes=args.shared_min_escape_celltypes,
        intermediate_min_escape_celltypes=args.intermediate_min_escape_celltypes,
        intermediate_max_escape_celltypes=args.intermediate_max_escape_celltypes,
        selective_max_escape_celltypes=args.selective_max_escape_celltypes,
        require_non_escape_for_modulated=args.require_non_escape_for_modulated,
        skip_03=args.skip_03,
    )


if __name__ == "__main__":
    main()
