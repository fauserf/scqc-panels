"""Plotting functions for QC diagnostics."""

from __future__ import annotations

import matplotlib.pyplot as plt
from anndata import AnnData

__all__ = ["qc_panel"]

# adata.obs column -> display label. Only columns present are plotted, so
# this works whether metrics came from scqc_panels.pp, scanpy, or both.
_DEFAULT_METRICS: dict[str, str] = {
    "n_genes_by_counts": "Genes per cell",
    "total_counts": "Total counts per cell",
    "pct_counts_mt": "Mitochondrial %",
    "pct_counts_ribo": "Ribosomal (cytosol) %",
    "pct_counts_malat1": "MALAT1 %",
    "doublet_score": "Doublet score",
}


def qc_panel(
    adata: AnnData,
    *,
    metrics: list[str] | None = None,
    ncols: int = 3,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Plot a standardized grid of QC diagnostic violin plots.

    Parameters
    ----------
    adata
        Annotated data matrix with QC metrics already computed in
        ``adata.obs`` (e.g. via :func:`scqc_panels.pp.compute_qc_metrics`
        and/or :func:`scanpy.pp.calculate_qc_metrics`).
    metrics
        Explicit list of ``adata.obs`` column names to plot, in order.
        Defaults to all recognized QC columns that are present.
    ncols
        Number of subplot columns in the grid.
    figsize
        Overall figure size. Defaults to a size scaled to the panel count.

    Returns
    -------
    The assembled :class:`matplotlib.figure.Figure`.
    """
    if metrics is None:
        metrics = [m for m in _DEFAULT_METRICS if m in adata.obs.columns]
    else:
        missing = [m for m in metrics if m not in adata.obs.columns]
        if missing:
            msg = f"Requested metrics not found in adata.obs: {missing!r}"
            raise ValueError(msg)

    if not metrics:
        msg = (
            "No recognized QC metrics found in adata.obs. Run "
            "scqc_panels.pp.compute_qc_metrics(adata) and/or "
            "scanpy.pp.calculate_qc_metrics(adata) first."
        )
        raise ValueError(msg)

    nrows = -(-len(metrics) // ncols)  # ceil division
    if figsize is None:
        figsize = (4 * ncols, 3.5 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    flat_axes = axes.ravel()

    for ax, metric in zip(flat_axes, metrics, strict=False):
        title = _DEFAULT_METRICS.get(metric, metric)
        ax.violinplot(adata.obs[metric].to_numpy(), showmedians=True)
        ax.set_xticks([])
        ax.set_title(title)
        ax.set_ylabel(title)

    for ax in flat_axes[len(metrics) :]:
        ax.axis("off")

    fig.suptitle("QC diagnostics", fontsize=14)
    fig.tight_layout()
    return fig
