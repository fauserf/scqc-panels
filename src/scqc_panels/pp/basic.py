"""Preprocessing / QC metric functions."""

from __future__ import annotations

import re
from typing import cast

import numpy as np
from anndata import AnnData

__all__ = ["add_malat1_fraction", "add_ribosomal_score", "compute_qc_metrics"]

# Matches ribosomal protein genes, e.g. RPS6, RPL10 (small/large subunit).
# Computed dynamically against adata.var_names rather than a hardcoded gene
# list, so it works regardless of species or exact naming convention.
_RIBO_PATTERN = re.compile(r"^RP[SL]\d", flags=re.IGNORECASE)


def add_malat1_fraction(
    adata: AnnData,
    *,
    malat1_names: tuple[str, ...] = ("MALAT1",),
    layer: str | None = None,
    inplace: bool = True,
) -> np.ndarray | None:
    """Compute the fraction of counts attributable to MALAT1 per cell.

    A high MALAT1 fraction is commonly used as a QC marker of nuclear or
    ambient RNA contamination and cell damage in single-cell RNA-seq.

    Parameters
    ----------
    adata
        Annotated data matrix. Genes are expected in ``adata.var_names``.
    malat1_names
        Gene symbol(s) to treat as MALAT1 (multiple names cover naming
        variants/aliases across references).
    layer
        Layer to use for counts. Defaults to ``adata.X``.
    inplace
        If ``True`` (default), store the result in
        ``adata.obs["pct_counts_malat1"]`` and return ``None``. If
        ``False``, return the array instead of modifying `adata`.

    Returns
    -------
    ``None`` if ``inplace=True``, otherwise a 1D array of per-cell
    percentages.
    """
    matrix = adata.X if layer is None else adata.layers[layer]
    if matrix is None:
        msg = "adata.X is None; pass `layer` pointing to a valid matrix."
        raise ValueError(msg)
    matrix = cast("np.ndarray", matrix)

    is_malat1 = adata.var_names.isin(malat1_names)
    if not is_malat1.any():
        msg = f"None of {malat1_names!r} found in adata.var_names."
        raise ValueError(msg)

    total_counts = np.asarray(matrix.sum(axis=1)).ravel()
    malat1_counts = np.asarray(matrix[:, is_malat1].sum(axis=1)).ravel()

    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(total_counts > 0, malat1_counts / total_counts * 100, 0.0)

    if inplace:
        adata.obs["pct_counts_malat1"] = pct
        return None
    return pct


def add_ribosomal_score(
    adata: AnnData,
    *,
    layer: str | None = None,
    inplace: bool = True,
) -> np.ndarray | None:
    r"""Compute the fraction of counts from cytosolic ribosomal protein genes.

    Genes are matched dynamically against ``adata.var_names`` using the
    pattern ``^RP[SL]\\d`` (e.g. ``RPS6``, ``RPL10``), covering both the
    small (RPS) and large (RPL) ribosomal subunit protein-coding genes.

    Parameters
    ----------
    adata
        Annotated data matrix.
    layer
        Layer to use for counts. Defaults to ``adata.X``.
    inplace
        If ``True`` (default), store the result in
        ``adata.obs["pct_counts_ribo"]`` (and mark matched genes in
        ``adata.var["ribo"]``) and return ``None``. If ``False``, return
        the array instead.

    Returns
    -------
    ``None`` if ``inplace=True``, otherwise a 1D array of per-cell
    percentages.
    """
    matrix = adata.X if layer is None else adata.layers[layer]
    if matrix is None:
        msg = "adata.X is None; pass `layer` pointing to a valid matrix."
        raise ValueError(msg)
    matrix = cast("np.ndarray", matrix)

    is_ribo = np.array([bool(_RIBO_PATTERN.match(name)) for name in adata.var_names])

    total_counts = np.asarray(matrix.sum(axis=1)).ravel()
    ribo_counts = np.asarray(matrix[:, is_ribo].sum(axis=1)).ravel() if is_ribo.any() else np.zeros_like(total_counts)

    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(total_counts > 0, ribo_counts / total_counts * 100, 0.0)

    if inplace:
        adata.var["ribo"] = is_ribo
        adata.obs["pct_counts_ribo"] = pct
        return None
    return pct


def compute_qc_metrics(
    adata: AnnData,
    *,
    malat1_names: tuple[str, ...] = ("MALAT1",),
    layer: str | None = None,
) -> None:
    """Compute all QC diagnostic metrics used by :func:`scqc_panels.pl.qc_panel`.

    Convenience wrapper running :func:`add_malat1_fraction` and
    :func:`add_ribosomal_score` in one call.

    Parameters
    ----------
    adata
        Annotated data matrix.
    malat1_names
        Gene symbol(s) to treat as MALAT1.
    layer
        Layer to use for counts. Defaults to ``adata.X``.
    """
    add_malat1_fraction(adata, malat1_names=malat1_names, layer=layer)
    add_ribosomal_score(adata, layer=layer)
