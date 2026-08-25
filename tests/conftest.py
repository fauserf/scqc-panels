import anndata as ad
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def adata() -> ad.AnnData:
    """A small synthetic AnnData with known counts for testing QC metrics."""
    var_names = ["MALAT1", "RPS6", "RPL10", "GENE_A", "GENE_B"]
    counts = np.array(
        [
            [10, 5, 5, 30, 50],  # total=100; malat1=10%, ribo=10%
            [0, 0, 0, 50, 50],  # total=100; malat1=0%,  ribo=0%
            [20, 20, 20, 20, 20],  # total=100; malat1=20%, ribo=40%
        ],
        dtype=np.float32,
    )
    return ad.AnnData(X=counts, var=pd.DataFrame(index=var_names))
