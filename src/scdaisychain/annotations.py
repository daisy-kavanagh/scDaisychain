from __future__ import annotations

from contextlib import contextmanager
from importlib import resources
from pathlib import Path
from typing import Iterator

import pandas as pd


ANNOTATION_FILES = {
    "human": "human_chrX_PAR_escape_annotation_for_xi_pipeline.tsv",
    "mouse": "mouse_chrX_PAR_escape_annotation_for_xi_pipeline.tsv",
}


def available_annotations() -> list[str]:
    """Return available built-in annotation names."""
    return sorted(ANNOTATION_FILES)


def _normalise_species(species: str) -> str:
    species = str(species).strip().lower()
    aliases = {
        "homo_sapiens": "human",
        "hsapiens": "human",
        "hg38": "human",
        "grch38": "human",
        "mus_musculus": "mouse",
        "mmusculus": "mouse",
        "mm39": "mouse",
        "grcm39": "mouse",
    }
    species = aliases.get(species, species)

    if species not in ANNOTATION_FILES:
        raise ValueError(
            f"Unknown annotation {species!r}. "
            f"Available annotations: {', '.join(available_annotations())}"
        )

    return species


def annotation_resource(species: str):
    """
    Return the importlib.resources object for a built-in annotation.

    This is useful internally. For user-facing filesystem paths, prefer
    annotation_path(species) as a context manager.
    """
    species = _normalise_species(species)
    filename = ANNOTATION_FILES[species]

    return resources.files("scdaisychain").joinpath(
        "resources",
        "annotations",
        filename,
    )


@contextmanager
def annotation_path(species: str) -> Iterator[Path]:
    """
    Yield a filesystem path to a packaged annotation file.

    Use as:

        with annotation_path("human") as path:
            ...

    The context-manager form is robust even if the package is installed from
    a wheel/zip-like source.
    """
    resource = annotation_resource(species)

    with resources.as_file(resource) as path:
        yield path


def load_annotation(species: str) -> pd.DataFrame:
    """Load a built-in annotation as a pandas DataFrame."""
    with annotation_path(species) as path:
        return pd.read_csv(path, sep="\t", dtype=str)
