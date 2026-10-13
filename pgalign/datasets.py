"""Toy datasets bundled with PGAlign for tutorials and tests."""

from importlib import resources

import numpy as np
import pandas as pd
import scanpy as sc


def load_mouse_pancreas_toy():
    """Return two AnnData objects (unknown, reference) from a downsampled mouse-pancreas MALDI slice.

    The data is a 120-spectrum, 601-bin slice of the public mouse pancreas dataset
    (Zenodo `10.5281/zenodo.3607915`, Prade et al., 2020, CC BY 4.0), drawn from the
    21 x 21 pixel window centred at (1636, 863) and the m/z range [100, 250] on a
    0.25 m/z grid. The "unknown" companion is the same spectra shifted two bins higher
    in m/z, with multiplicative lognormal noise (sigma 0.05), so PGAlign should recover
    a shift of +2 bins. ``scripts/make_toy_data.py`` in the repository rebuilds it.

    Returns
    -------
    unk : anndata.AnnData
        The shifted "unknown" spectrum.
    ref : anndata.AnnData
        The unshifted reference spectrum.
    """
    with resources.files(__package__).joinpath("data/mouse_pancreas_toy.npz").open("rb") as f:
        npz = np.load(f)
        X_ref = npz["X_ref"]
        X_unk = npz["X_unk"]
        mz = npz["mz"]
        x = npz["x"]
        y = npz["y"]

    var = pd.DataFrame({"m/z": mz.astype(float)}, index=[str(v) for v in mz])
    obs = pd.DataFrame({"x": x, "y": y}, index=[f"s{i}" for i in range(len(x))])
    ref = sc.AnnData(X=X_ref, var=var, obs=obs)
    unk = sc.AnnData(X=X_unk, var=var.copy(), obs=obs.copy())
    return unk, ref
