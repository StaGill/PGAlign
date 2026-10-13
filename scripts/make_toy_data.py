"""Rebuild pgalign/data/mouse_pancreas_toy.npz from the gridded mouse pancreas data.

Source: MousePancreate_small.h5ad, the 0.25 m/z gridded version of
Zenodo 10.5281/zenodo.3607915 (Prade et al., 2020, CC BY 4.0) produced by
CodeInPaper/MousePancreate.ipynb.

Usage: python scripts/make_toy_data.py path/to/MousePancreate_small.h5ad
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

CENTER = (1636, 863)
HALF_WIDTH = 10
N_SPECTRA = 120
MZ_RANGE = (100.0, 250.0)
SHIFT_BINS = 2
NOISE_SIGMA = 0.05
SEED = 0

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "pgalign" / "data" / "mouse_pancreas_toy.npz"


def main():
    parser = argparse.ArgumentParser(description="Rebuild the bundled mouse pancreas toy data.")
    parser.add_argument("source", help="path to MousePancreate_small.h5ad")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="output .npz path")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.source, backed="r")
    x = adata.obs["x"].to_numpy()
    y = adata.obs["y"].to_numpy()
    window = np.flatnonzero((np.abs(x - CENTER[0]) <= HALF_WIDTH) & (np.abs(y - CENTER[1]) <= HALF_WIDTH))
    rng = np.random.default_rng(SEED)
    rows = np.sort(rng.choice(window, size=N_SPECTRA, replace=False))

    mz = adata.var["m/z"].to_numpy(dtype=float)
    cols = np.flatnonzero((mz >= MZ_RANGE[0]) & (mz <= MZ_RANGE[1]))
    if cols[0] < SHIFT_BINS:
        raise ValueError("the source grid must extend SHIFT_BINS bins below the m/z range")

    X = adata[rows, :].to_memory().X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X_ref = X[:, cols].astype(np.float32)
    # Reading every spectrum SHIFT_BINS bins lower in m/z moves each peak SHIFT_BINS bins higher.
    noise = rng.lognormal(mean=0.0, sigma=NOISE_SIGMA, size=X_ref.shape)
    X_unk = (X[:, cols - SHIFT_BINS] * noise).astype(np.float32)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        X_ref=X_ref,
        X_unk=X_unk,
        mz=mz[cols],
        x=x[rows].astype(np.int64),
        y=y[rows].astype(np.int64),
    )
    adata.file.close()
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KB): {X_ref.shape[0]} spectra x {X_ref.shape[1]} bins")


if __name__ == "__main__":
    main()
