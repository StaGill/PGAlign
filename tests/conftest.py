"""Shared fixtures for GALAXY smoke tests."""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc


def _make_anndata(mz_axis, n_spectra, peak_centers, peak_widths, rng):
    """Build an AnnData object with synthetic Gaussian peaks."""
    n_bins = len(mz_axis)
    X = np.zeros((n_spectra, n_bins))
    for center, width in zip(peak_centers, peak_widths):
        # Gaussian peak centred at ``center`` with sigma ``width``.
        peak = np.exp(-0.5 * ((mz_axis - center) / width) ** 2)
        # Per-spectrum amplitude jitter so the peak height varies a little.
        amps = 1.0 + 0.05 * rng.standard_normal(n_spectra)
        X += np.outer(amps, peak)
    # Mild background noise so the derivative-based peak caller has something to threshold.
    X += 0.01 * rng.standard_normal(X.shape)
    var = pd.DataFrame({"m/z": mz_axis}, index=[str(v) for v in mz_axis])
    obs = pd.DataFrame(index=[f"s{i}" for i in range(n_spectra)])
    return sc.AnnData(X=X, var=var, obs=obs)


@pytest.fixture
def tiny_anndata():
    """Return two AnnData objects with known peaks, the unknown shifted by 1 bin."""
    rng = np.random.default_rng(seed=0)
    increment = 0.5
    mz_axis = np.arange(100.0, 200.0 + increment, increment)
    peak_centers = [120.0, 140.0, 160.0, 180.0]
    peak_widths = [0.8] * len(peak_centers)
    ref = _make_anndata(mz_axis, n_spectra=20, peak_centers=peak_centers,
                        peak_widths=peak_widths, rng=rng)
    # Unknown spectrum: same peaks but shifted by ``increment`` m/z units.
    unk = _make_anndata(mz_axis, n_spectra=20,
                        peak_centers=[c + increment for c in peak_centers],
                        peak_widths=peak_widths, rng=rng)
    return unk, ref
