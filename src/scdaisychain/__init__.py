"""scDaisychain XCI pipeline package."""

from .anndata_io import (
    add_x_layers_and_metrics,
    add_haplotype_sums,
    load_10x_batches_with_optional_layers,
)
from .xi_downstream import run_xi_downstream_from_adata
from .annotations import (
    annotation_path,
    available_annotations,
    load_annotation,
)
__version__ = "0.1.0"

__all__ = [
    "add_x_layers_and_metrics",
    "add_haplotype_sums",
    "load_10x_batches_with_optional_layers",
    "run_xi_downstream_from_adata",
]