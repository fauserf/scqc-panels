import anndata as ad
import numpy as np
import pytest

import scqc_panels


def test_add_malat1_fraction(adata):
    scqc_panels.pp.add_malat1_fraction(adata)

    assert "pct_counts_malat1" in adata.obs.columns
    np.testing.assert_allclose(adata.obs["pct_counts_malat1"].to_numpy(), [10.0, 0.0, 20.0])


def test_add_malat1_fraction_missing_gene(adata):
    without_malat1 = adata[:, adata.var_names != "MALAT1"].copy()
    with pytest.raises(ValueError, match="MALAT1"):
        scqc_panels.pp.add_malat1_fraction(without_malat1)


def test_add_ribosomal_score(adata):
    scqc_panels.pp.add_ribosomal_score(adata)

    assert "pct_counts_ribo" in adata.obs.columns
    assert adata.var["ribo"].tolist() == [False, True, True, False, False]
    np.testing.assert_allclose(adata.obs["pct_counts_ribo"].to_numpy(), [10.0, 0.0, 40.0])


def test_compute_qc_metrics(adata):
    scqc_panels.pp.compute_qc_metrics(adata)

    assert "pct_counts_malat1" in adata.obs.columns
    assert "pct_counts_ribo" in adata.obs.columns


def test_qc_panel_returns_figure(adata):
    scqc_panels.pp.compute_qc_metrics(adata)
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()

    fig = scqc_panels.pl.qc_panel(adata)
    assert fig is not None


def test_qc_panel_raises_without_metrics():
    empty_adata = ad.AnnData(X=np.zeros((3, 2), dtype=np.float32))
    with pytest.raises(ValueError, match="No recognized QC metrics"):
        scqc_panels.pl.qc_panel(empty_adata)
