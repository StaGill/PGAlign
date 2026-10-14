"""Readers that turn MSI exports into AnnData objects for :func:`pgalign.align`, and .h5ad I/O."""

import math

import numpy as np
import pandas as pd
import scanpy as sc

from .util import gridding


def read_imzml(path, increment=None, mz_range=None, epsilon=0.025):
    """Read an imzML file into an AnnData that ``pg.align`` accepts.

    In continuous mode every spectrum shares one m/z axis, which is used as it is.
    In processed mode each spectrum carries its own m/z values, so ``increment`` is
    required and every spectrum is put on a common grid with ``gridding``
    (Nadaraya-Watson smoothing with bandwidth ``epsilon``), as the manuscript's
    simulation study did with ``increment=0.25``. Passing ``increment`` grids the
    spectra in either mode. Requires pyimzML: ``pip install "pgalign[imzml]"``.

    :param path: the ``.imzML`` file; its ``.ibd`` file must sit beside it.
    :param increment: grid spacing in m/z. Required when the spectra do not share
        one m/z axis.
    :param mz_range: ``(start, end)`` of the grid, end excluded. Defaults to one
        below the floor of the smallest m/z in the file and one above the ceiling
        of the largest, as in ``gridding``.
    :param epsilon: bandwidth of the gridding kernel. Defaults to 0.025.
    :return: an ``anndata.AnnData`` with the m/z values in ``var["m/z"]`` and the
        pixel coordinates in ``obs`` (``x``, ``y``, ``z``).
    """
    try:
        from pyimzml.ImzMLParser import ImzMLParser
    except ImportError as err:
        raise ImportError('read_imzml needs pyimzML: pip install "pgalign[imzml]"') from err
    parser = ImzMLParser(path)
    spectra = [parser.getspectrum(i) for i in range(len(parser.coordinates))]
    if increment is None:
        mz = np.asarray(spectra[0][0], dtype=float)
        if not all(np.array_equal(mzs, mz) for mzs, _ in spectra):
            raise ValueError(
                f"the spectra in {path} do not share one m/z axis (processed imzML); "
                "pass increment= to grid them, for example increment=0.25 as in the "
                "manuscript's simulation study"
            )
        X = np.vstack([np.asarray(intensities, dtype=float) for _, intensities in spectra])
    else:
        observed = [mzs for mzs, _ in spectra if len(mzs)]
        lo, hi = min(m.min() for m in observed), max(m.max() for m in observed)
        if mz_range is None:
            mz_range = (math.floor(lo) - 1, math.ceil(hi) + 1)
        elif mz_range[0] > hi or mz_range[1] <= lo:
            raise ValueError(
                f"no m/z value in {path} falls inside mz_range {tuple(mz_range)}; "
                f"the file covers {lo:g} to {hi:g}"
            )
        grids = [gridding(np.asarray(mzs, dtype=float), np.asarray(intensities, dtype=float),
                          start=mz_range[0], end=mz_range[1], increment=increment, epsilon=epsilon)
                 for mzs, intensities in spectra]
        mz = grids[0][0]
        X = np.vstack([grid[1] for grid in grids])
    var = pd.DataFrame({"m/z": mz}, index=[str(v) for v in mz])
    obs = pd.DataFrame(np.asarray(parser.coordinates)[:, :3], columns=["x", "y", "z"],
                       index=[f"s{i}" for i in range(len(spectra))])
    return sc.AnnData(X=X, var=var, obs=obs)


def write_h5ad(adata, path):
    """Write an AnnData to ``.h5ad``, storing ``var["m/z"]`` as ``var["mz"]``.

    anndata 0.13 refuses to write a column whose name contains a slash, so the m/z
    values are stored under ``mz``; :func:`read_h5ad` restores them. The object
    passed in is not modified, and its data are written without being copied.

    :param adata: the ``anndata.AnnData`` to write, for example ``res.aligned``.
    :param path: the ``.h5ad`` file to write.
    """
    if "m/z" in adata.var.columns and "mz" in adata.var.columns:
        raise ValueError("var has both an 'm/z' and an 'mz' column; drop one before writing")
    out = sc.AnnData(
        X=adata.X, obs=adata.obs, var=adata.var.rename(columns={"m/z": "mz"}),
        uns=dict(adata.uns), obsm=dict(adata.obsm), varm=dict(adata.varm),
        layers=dict(adata.layers), obsp=dict(adata.obsp), varp=dict(adata.varp),
    )
    out.write_h5ad(path)


def read_h5ad(path):
    """Read an ``.h5ad`` file written by :func:`write_h5ad` or by older anndata.

    Restores ``var["mz"]`` to ``var["m/z"]`` when the file has no ``m/z`` column;
    files that store ``m/z`` itself are read as they are.

    :param path: the ``.h5ad`` file to read.
    :return: an ``anndata.AnnData`` with the m/z values in ``var["m/z"]``.
    """
    adata = sc.read_h5ad(path)
    if "m/z" not in adata.var.columns and "mz" in adata.var.columns:
        adata.var = adata.var.rename(columns={"mz": "m/z"})
    return adata


def read_scils_csv(spectra_csv, spots_csv=None, sep=";"):
    """Read a SCiLS Lab CSV export into an AnnData that ``pg.align`` accepts.

    The spectra file holds one row per spot: the spot name in the first column,
    then one intensity column per m/z value, with the m/z values as the column
    headers. The optional spots file (``Region Spots.csv``) holds one row per spot,
    in the same order, with its index and x and y coordinates; its columns are
    copied into ``obs`` unchanged.

    :param spectra_csv: path of the spectra export (``... - All Spectra.csv``).
    :param spots_csv: path of the spot coordinates (``... - Region Spots.csv``).
    :param sep: the field separator. Defaults to ``";"``, as SCiLS writes it.
    :return: an ``anndata.AnnData`` with the m/z values in ``var["m/z"]``.
    """
    spectra = pd.read_csv(spectra_csv, sep=sep)
    labels = list(spectra.columns[1:])
    var = pd.DataFrame({"m/z": np.asarray(labels, dtype=float)}, index=labels)
    obs = pd.DataFrame(index=spectra.iloc[:, 0].astype(str).to_numpy())
    if spots_csv is not None:
        spots = pd.read_csv(spots_csv, sep=sep)
        if len(spots) != len(spectra):
            raise ValueError(
                f"{spots_csv} has {len(spots)} spots but {spectra_csv} has {len(spectra)} spectra"
            )
        for column in spots.columns:
            obs[column] = spots[column].to_numpy()
    X = spectra.iloc[:, 1:].to_numpy(dtype=float)
    return sc.AnnData(X=X, var=var, obs=obs)
