"""Input checks of pgalign.align."""

import numpy as np
import pytest
import scipy.sparse as sp

import pgalign as pg


@pytest.fixture(scope="module")
def toy():
    return pg.datasets.load_mouse_pancreas_toy()


def _renamed(adata):
    adata = adata.copy()
    adata.var = adata.var.rename(columns={"m/z": "mz"})
    return adata


def _descending(adata):
    return adata[:, ::-1].copy()


def _with_nan(adata):
    adata = adata.copy()
    X = np.asarray(adata.X).copy()
    X[0, 5] = np.nan
    adata.X = X
    return adata


def _moved(adata):
    adata = adata.copy()
    adata.var["m/z"] = adata.var["m/z"] + 1000
    return adata


@pytest.mark.parametrize("damage, message", [
    (_renamed, "has no 'm/z' column"),
    (_descending, "must be strictly increasing"),
    (_with_nan, "1 missing or infinite intensities"),
    (_moved, "m/z ranges do not overlap"),
])
def test_unusable_input_raises_with_the_fix(toy, damage, message):
    unk, ref = toy
    with pytest.raises(ValueError, match=message):
        pg.align(damage(unk), ref)


def test_sparse_input_gives_the_dense_result(toy):
    unk, ref = toy
    sparse_unk, sparse_ref = unk.copy(), ref.copy()
    sparse_unk.X = sp.csr_matrix(sparse_unk.X)
    sparse_ref.X = sp.csr_matrix(sparse_ref.X)
    res = pg.align(sparse_unk, sparse_ref)
    np.testing.assert_array_equal(np.asarray(res.aligned.X), np.asarray(pg.align(unk, ref).aligned.X))
    assert sp.issparse(sparse_unk.X)
