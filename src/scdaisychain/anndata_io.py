from __future__ import annotations

"""
AnnData / 10x loading helpers for scDaisychain.

This module provides convenience functions for loading one or more 10x matrices
and optionally attaching scDaisychain Xa/Xi count layers and haplotype-level QC
metrics.

Main public functions
---------------------
- add_x_layers_and_metrics
- add_haplotype_sums
- load_10x_batches_with_optional_layers
"""

from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


PathLike = Union[str, Path]
MatrixSpec = Union[
    PathLike,
    Tuple[PathLike, str],
    Sequence[Union[PathLike, Tuple[PathLike, str]]],
]


def _read_optional_table(
    path: Optional[PathLike],
    *,
    index_col: Optional[int] = 0,
    sep: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Read an optional CSV/TSV-like table."""
    if path is None:
        return None

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")

    if sep is None:
        # Let pandas infer comma vs tab. This is slower but robust for small
        # side-car tables like Xa/Xi matrices and haplotype sums.
        return pd.read_csv(p, index_col=index_col, sep=None, engine="python")

    return pd.read_csv(p, index_col=index_col, sep=sep)


def _as_csr_float(X) -> sparse.csr_matrix:
    """Return X as float CSR sparse matrix."""
    if sparse.issparse(X):
        return X.tocsr().astype(float, copy=False)
    return sparse.csr_matrix(np.asarray(X, dtype=float))


def _as_dense_float(X) -> np.ndarray:
    """Return X as dense float ndarray."""
    if sparse.issparse(X):
        return X.toarray().astype(float, copy=False)
    return np.asarray(X, dtype=float)


def _replace_nan_inplace(X):
    """Replace NaNs in sparse/dense matrix-like object and return the object."""
    if sparse.issparse(X):
        X = X.tocsr(copy=False).astype(float, copy=False)
        if X.data.size:
            X.data[np.isnan(X.data)] = 0.0
            X.eliminate_zeros()
        return X

    X = np.asarray(X, dtype=float)
    np.nan_to_num(X, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    return X


def _sum_axis(X, axis: int) -> np.ndarray:
    """Sum dense/sparse matrix along an axis and return a 1D ndarray."""
    if sparse.issparse(X):
        return np.asarray(X.sum(axis=axis)).ravel()
    return np.asarray(X.sum(axis=axis)).ravel()


def _row_scale(X, scale: np.ndarray):
    """Multiply each row of X by scale, preserving sparse when possible."""
    scale = np.asarray(scale, dtype=float)
    scale = np.nan_to_num(scale, nan=0.0, posinf=0.0, neginf=0.0)

    if sparse.issparse(X):
        return X.multiply(scale[:, None]).tocsr()

    return np.asarray(X, dtype=float) * scale[:, None]


def _add_suffix_if_needed(values: Sequence[str], suffix: Optional[str]) -> List[str]:
    """Append suffix to values unless already present. If suffix is None, do nothing."""
    if suffix is None or suffix == "":
        return [str(v) for v in values]

    out = []
    for v in values:
        s = str(v)
        out.append(s if s.endswith(suffix) else f"{s}{suffix}")
    return out


def _align_gene_by_cell_table_to_adata(
    table: pd.DataFrame,
    adata: ad.AnnData,
    *,
    layer_name: str,
    barcode_suffix: Optional[str] = None,
    fill_missing: float = 0.0,
    make_dense: bool = False,
    warn: bool = True,
):
    """
    Align a gene x cell table to AnnData cells x genes.

    The input table is expected to have:
      - rows = gene names matching adata.var_names
      - columns = cell barcodes matching adata.obs_names

    Returns a cells x genes matrix suitable for adata.layers[layer_name].
    """
    df = table.copy()
    df.index = df.index.astype(str)
    df.columns = _add_suffix_if_needed(df.columns.astype(str).tolist(), barcode_suffix)

    obs_names = pd.Index(adata.obs_names.astype(str))
    var_names = pd.Index(adata.var_names.astype(str))

    n_gene_matches = int(var_names.isin(df.index).sum())
    n_cell_matches = int(obs_names.isin(df.columns).sum())

    if warn:
        if n_gene_matches == 0:
            warnings.warn(
                f"No genes in {layer_name} table matched adata.var_names. "
                "Check whether the table uses gene symbols vs Ensembl IDs.",
                stacklevel=2,
            )
        if n_cell_matches == 0:
            warnings.warn(
                f"No cells in {layer_name} table matched adata.obs_names. "
                "Check barcode suffixes such as '-1'.",
                stacklevel=2,
            )

    aligned = df.reindex(index=var_names, columns=obs_names, fill_value=fill_missing)
    aligned = aligned.apply(pd.to_numeric, errors="coerce").fillna(fill_missing)

    # Transpose from genes x cells to cells x genes.
    values = aligned.to_numpy(dtype=float).T
    return values if make_dense else sparse.csr_matrix(values)


def add_x_layers_and_metrics(
    adata: ad.AnnData,
    xa_csv: Optional[PathLike] = None,
    xi_csv: Optional[PathLike] = None,
    *,
    xa_layer: str = "Xa",
    xi_layer: str = "Xi",
    xai_layer: str = "Xai",
    raw_layer: str = "raw_counts",
    make_dense_layers: bool = False,
    normalize_target_sum: float = 10_000.0,
    barcode_suffix: Optional[str] = None,
    copy: bool = False,
) -> ad.AnnData:
    """
    Add Xa/Xi layers and derived chrX allelic metrics to an AnnData object.

    Parameters
    ----------
    adata
        AnnData object with cells in rows and genes in columns.
    xa_csv / xi_csv
        Optional gene x cell CSV/TSV tables. Rows should match adata.var_names
        and columns should match adata.obs_names. If your CSV barcodes lack a
        suffix present in adata.obs_names, use barcode_suffix="-1".
    xa_layer / xi_layer
        Names of the layers to create/use.
    xai_layer
        Name for combined Xa+Xi layer.
    raw_layer
        Name of layer storing a copy of the input adata.X. If absent, it is
        created as adata.X.copy().
    make_dense_layers
        If False, Xa/Xi/Xai and normalized layers are stored as CSR sparse
        matrices. This is safer for large single-cell data. Set True only for
        small matrices or when dense arrays are explicitly desired.
    normalize_target_sum
        Per-cell scale target for Xa/Xi normalized layers. Scaling denominator
        is raw_total from raw_layer.
    barcode_suffix
        Optional suffix to append to Xa/Xi CSV columns before alignment, e.g. "-1".
        Suffix is only appended if not already present.
    copy
        If True, returns a modified copy. Otherwise modifies adata in place.

    Adds
    ----
    Layers:
      - raw_layer
      - xa_layer
      - xi_layer
      - xai_layer
      - f"{xa_layer}_Normalized"
      - f"{xi_layer}_Normalized"
      - f"{xai_layer}_Normalized"

    obs columns:
      - raw_total
      - Xai_total
      - allelic_imbalance_chrX
    """
    if copy:
        adata = adata.copy()

    if raw_layer not in adata.layers:
        adata.layers[raw_layer] = adata.X.copy()

    xa_df = _read_optional_table(xa_csv, index_col=0)
    xi_df = _read_optional_table(xi_csv, index_col=0)

    if xa_df is not None:
        adata.layers[xa_layer] = _align_gene_by_cell_table_to_adata(
            xa_df,
            adata,
            layer_name=xa_layer,
            barcode_suffix=barcode_suffix,
            fill_missing=0.0,
            make_dense=make_dense_layers,
        )

    if xi_df is not None:
        adata.layers[xi_layer] = _align_gene_by_cell_table_to_adata(
            xi_df,
            adata,
            layer_name=xi_layer,
            barcode_suffix=barcode_suffix,
            fill_missing=0.0,
            make_dense=make_dense_layers,
        )

    if xa_layer not in adata.layers or xi_layer not in adata.layers:
        return adata

    # Clean NaNs and standardise storage.
    if make_dense_layers:
        Xa = _as_dense_float(_replace_nan_inplace(adata.layers[xa_layer]))
        Xi = _as_dense_float(_replace_nan_inplace(adata.layers[xi_layer]))
    else:
        Xa = _as_csr_float(_replace_nan_inplace(adata.layers[xa_layer]))
        Xi = _as_csr_float(_replace_nan_inplace(adata.layers[xi_layer]))

    adata.layers[xa_layer] = Xa
    adata.layers[xi_layer] = Xi
    adata.layers[xai_layer] = Xa + Xi

    raw_counts = adata.layers[raw_layer]
    adata.obs["raw_total"] = _sum_axis(raw_counts, axis=1)

    denom = adata.obs["raw_total"].to_numpy(dtype=float)
    scale = np.divide(
        float(normalize_target_sum),
        denom,
        out=np.zeros_like(denom, dtype=float),
        where=denom > 0,
    )

    xa_norm_layer = f"{xa_layer}_Normalized"
    xi_norm_layer = f"{xi_layer}_Normalized"
    xai_norm_layer = f"{xai_layer}_Normalized"

    adata.layers[xa_norm_layer] = _row_scale(Xa, scale)
    adata.layers[xi_norm_layer] = _row_scale(Xi, scale)
    adata.layers[xai_norm_layer] = adata.layers[xa_norm_layer] + adata.layers[xi_norm_layer]

    adata.obs["Xai_total"] = _sum_axis(adata.layers[xai_norm_layer], axis=1)

    xa_sum = _sum_axis(adata.layers[xa_norm_layer], axis=1)
    xi_sum = _sum_axis(adata.layers[xi_norm_layer], axis=1)
    denom2 = xa_sum + xi_sum

    with np.errstate(divide="ignore", invalid="ignore"):
        frac_xa = np.divide(
            xa_sum,
            denom2,
            out=np.full_like(denom2, np.nan, dtype=float),
            where=denom2 > 0,
        )
        adata.obs["allelic_imbalance_chrX"] = 2 * np.abs(0.5 - frac_xa)

    return adata


def add_haplotype_sums(
    adata: ad.AnnData,
    haplotype_sums_csv: Optional[PathLike] = None,
    *,
    cell_barcode_col: str = "cell_barcode",
    hap1_col: str = "Haplotype_1_Sum",
    hap2_col: str = "Haplotype_2_Sum",
    active_x_col: str = "active_X",
    suffix_barcodes: Optional[str] = "-1",
    copy: bool = False,
) -> ad.AnnData:
    """
    Add haplotype-level sums and QC metrics to adata.obs.

    Parameters
    ----------
    haplotype_sums_csv
        Table containing one row per cell barcode and at least hap1_col/hap2_col.
    suffix_barcodes
        Optional suffix to append to CSV barcodes before matching adata.obs_names,
        e.g. "-1". Suffix is only appended if not already present. Set None to
        disable suffix modification.

    Adds
    ----
    obs columns:
      - active_X, if present in the CSV
      - allelic_imbalance_HC
      - haplotype_sum
    """
    if copy:
        adata = adata.copy()

    df = _read_optional_table(haplotype_sums_csv, index_col=None)
    if df is None:
        return adata

    required = {cell_barcode_col, hap1_col, hap2_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"haplotype_sums_csv is missing required columns: {sorted(missing)}")

    df = df.copy()
    df[hap1_col] = pd.to_numeric(df[hap1_col], errors="coerce")
    df[hap2_col] = pd.to_numeric(df[hap2_col], errors="coerce")

    denom = df[hap1_col].to_numpy(dtype=float) + df[hap2_col].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.divide(
            df[hap1_col].to_numpy(dtype=float),
            denom,
            out=np.full_like(denom, np.nan, dtype=float),
            where=denom > 0,
        )
        df["allelic_imbalance"] = 2 * np.abs(0.5 - frac)

    df["haplotype_sum"] = denom
    df[cell_barcode_col] = _add_suffix_if_needed(df[cell_barcode_col].astype(str).tolist(), suffix_barcodes)

    if df[cell_barcode_col].duplicated().any():
        duplicated = df.loc[df[cell_barcode_col].duplicated(), cell_barcode_col].head(5).tolist()
        warnings.warn(
            f"Duplicate cell barcodes in haplotype_sums_csv after suffix handling. "
            f"Keeping first occurrence. Examples: {duplicated}",
            stacklevel=2,
        )
        df = df.drop_duplicates(subset=[cell_barcode_col], keep="first")

    df = df.set_index(cell_barcode_col)
    aligned = df.reindex(index=adata.obs_names)

    if active_x_col in aligned.columns:
        adata.obs["active_X"] = aligned[active_x_col].to_numpy()

    adata.obs["allelic_imbalance_HC"] = aligned["allelic_imbalance"].to_numpy()
    adata.obs["haplotype_sum"] = aligned["haplotype_sum"].to_numpy()

    n_matched = int(pd.Index(adata.obs_names).isin(df.index).sum())
    if n_matched == 0:
        warnings.warn(
            "No haplotype_sums_csv barcodes matched adata.obs_names. "
            "Check suffix_barcodes.",
            stacklevel=2,
        )

    return adata


def _normalise_matrix_specs(matrices: MatrixSpec) -> List[Tuple[Path, str]]:
    """Normalise supported input forms into [(matrix_dir, batch_name), ...]."""
    if isinstance(matrices, (str, Path)):
        return [(Path(matrices), "batch1")]

    matrices_list = list(matrices)  # type: ignore[arg-type]
    if not matrices_list:
        raise ValueError("You must provide at least one matrix directory.")

    batch_items: List[Tuple[Path, str]] = []
    auto_i = 1
    for item in matrices_list:
        if isinstance(item, (tuple, list)) and len(item) == 2:
            d, name = item
            batch_items.append((Path(d), str(name)))
        else:
            batch_items.append((Path(item), f"batch{auto_i}"))
            auto_i += 1

    return batch_items


def load_10x_batches_with_optional_layers(
    matrices: MatrixSpec,
    *,
    var_names: str = "gene_symbols",
    cache: bool = False,
    sample: Optional[Union[str, int, float]] = None,
    sample_by_batch: Optional[Mapping[str, Union[str, int, float]]] = None,
    xa_csv_by_batch: Optional[Mapping[str, PathLike]] = None,
    xi_csv_by_batch: Optional[Mapping[str, PathLike]] = None,
    haplotype_sums_csv: Optional[PathLike] = None,
    haplotype_sums_csv_by_batch: Optional[Mapping[str, PathLike]] = None,
    join: str = "outer",
    batch_obs_key: str = "batch",
    index_unique: Optional[str] = "-",
    normalize_target_sum: float = 10_000.0,
    make_dense_layers: bool = False,
    barcode_suffix: Optional[str] = None,
    haplotype_barcode_suffix: Optional[str] = "-1",
) -> ad.AnnData:
    """
    Load one or more 10x matrix directories and optionally attach Xa/Xi layers.

    Parameters
    ----------
    matrices
        One of:
          - "path/to/10x_dir"
          - ("path/to/10x_dir", "batch_name")
          - [("dir1", "name1"), ("dir2", "name2"), ...]
          - ["dir1", "dir2", ...] with auto batch names batch1, batch2, ...
    sample
        One value assigned to adata.obs["sample"] for every batch.
    sample_by_batch
        Optional mapping batch_name -> sample label. Overrides sample for those
        batches.
    xa_csv_by_batch / xi_csv_by_batch
        Mapping batch_name -> Xa/Xi CSV. CSVs should be gene x cell.
    haplotype_sums_csv
        One haplotype sums CSV applied to all batches.
    haplotype_sums_csv_by_batch
        Mapping batch_name -> haplotype sums CSV. Overrides shared
        haplotype_sums_csv for those batches.
    join
        Passed to anndata.concat.
    batch_obs_key
        obs column for batch labels.
    index_unique
        Passed to anndata.concat. Default "-" makes duplicate barcodes unique
        across batches, e.g. AAAC-1-batch1. Set None if obs_names are already
        globally unique.
    make_dense_layers
        If False, Xa/Xi layers are sparse CSR. Recommended for large data.
    barcode_suffix
        Optional suffix to append to Xa/Xi CSV cell columns before matching.
        Use "-1" if CSV barcodes lack the Cell Ranger suffix.
    haplotype_barcode_suffix
        Optional suffix to append to haplotype CSV barcodes before matching.
        Defaults to "-1" to preserve the behaviour of the original helper.

    Returns
    -------
    AnnData
        Merged AnnData object with optional Xa/Xi layers and QC metrics.
    """
    batch_items = _normalise_matrix_specs(matrices)

    xa_csv_by_batch = dict(xa_csv_by_batch or {})
    xi_csv_by_batch = dict(xi_csv_by_batch or {})
    haplotype_sums_csv_by_batch = dict(haplotype_sums_csv_by_batch or {})
    sample_by_batch = dict(sample_by_batch or {})

    adatas: List[ad.AnnData] = []
    keys: List[str] = []

    for matrix_dir, batch_name in batch_items:
        if not matrix_dir.exists():
            raise FileNotFoundError(f"10x matrix directory not found: {matrix_dir}")

        a = sc.read_10x_mtx(str(matrix_dir), var_names=var_names, cache=cache)

        batch_sample = sample_by_batch.get(batch_name, sample)
        if batch_sample is not None:
            a.obs["sample"] = batch_sample

        a.obs[batch_obs_key] = batch_name

        a = add_x_layers_and_metrics(
            a,
            xa_csv=xa_csv_by_batch.get(batch_name),
            xi_csv=xi_csv_by_batch.get(batch_name),
            make_dense_layers=make_dense_layers,
            normalize_target_sum=normalize_target_sum,
            barcode_suffix=barcode_suffix,
        )

        hs = haplotype_sums_csv_by_batch.get(batch_name, haplotype_sums_csv)
        a = add_haplotype_sums(
            a,
            haplotype_sums_csv=hs,
            suffix_barcodes=haplotype_barcode_suffix,
        )

        adatas.append(a)
        keys.append(batch_name)

    if len(adatas) == 1:
        out = adatas[0]
        out.obs[batch_obs_key] = keys[0]
        return out

    merged = ad.concat(
        adatas,
        join=join,
        label=batch_obs_key,
        keys=keys,
        index_unique=index_unique,
    )

    if sample is not None and "sample" not in merged.obs:
        merged.obs["sample"] = sample

    return merged


# Backward-compatible private aliases, in case notebooks already used these names.
_add_x_layers_and_metrics = add_x_layers_and_metrics
_add_optional_haplotype_sums = add_haplotype_sums


__all__ = [
    "add_x_layers_and_metrics",
    "add_haplotype_sums",
    "load_10x_batches_with_optional_layers",
]
