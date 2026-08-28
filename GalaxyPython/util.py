"""Utility helpers used by GALAXY's peak-calling and alignment routines."""

import logging
import math

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


def prefilter_cells(adata, min_counts=None, max_counts=None, min_genes=200, max_genes=None):
    if min_genes is None and min_counts is None and max_genes is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp = np.asarray([True] * adata.shape[0], dtype=bool)
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_cells(adata.X, min_genes=min_genes)[0]) if min_genes is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_cells(adata.X, max_genes=max_genes)[0]) if max_genes is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_cells(adata.X, min_counts=min_counts)[0]) if min_counts is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_cells(adata.X, max_counts=max_counts)[0]) if max_counts is not None else id_tmp
    adata._inplace_subset_obs(id_tmp)
    adata.raw = sc.pp.log1p(adata, copy=True)
    logger.debug("adata.raw.var_names.is_unique = %s", adata.raw.var_names.is_unique)


def prefilter_genes(adata, min_counts=None, max_counts=None, min_cells=10, max_cells=None):
    if min_cells is None and min_counts is None and max_cells is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp = np.asarray([True] * adata.shape[1], dtype=bool)
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, min_cells=min_cells)[0]) if min_cells is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, max_cells=max_cells)[0]) if max_cells is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, min_counts=min_counts)[0]) if min_counts is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, max_counts=max_counts)[0]) if max_counts is not None else id_tmp
    adata._inplace_subset_var(id_tmp)


def prefilter_specialgenes(adata, Gene1Pattern="ERCC", Gene2Pattern="MT-"):
    id_tmp1 = np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names], dtype=bool)
    id_tmp2 = np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names], dtype=bool)
    id_tmp = np.logical_and(id_tmp1, id_tmp2)
    adata._inplace_subset_var(id_tmp)


def save_clusterresults(cluster, filepath):
    rangetxt = np.zeros((len(cluster), 2))
    for i, points in enumerate(cluster):
        rangetxt[i, 0] = points[0]
        rangetxt[i, 1] = points[-1]
    np.savetxt("{filepath}.csv".format(filepath=filepath), rangetxt, delimiter=",")


def rotate(l, n):
    return l[n:] + l[:n]


def group_range(cluster):
    lows = np.array([i[0] for i in cluster])
    ups = np.array([i[-1] for i in cluster])
    return [lows, ups]


def in_range_lookup(mzvalues, cluster, expand=True):
    lows, ups = group_range(cluster)
    TFlist = [np.any((lows <= x) & (x <= ups)) for x in mzvalues]
    if expand:
        TFlistr = rotate(TFlist, 1)
        TFlistl = rotate(TFlist, -1)
        newTFlist = [any((mz1, mz2, mz3)) for (mz1, mz2, mz3) in zip(TFlist, TFlistr, TFlistl)]
    else:
        newTFlist = TFlist
    return newTFlist


def comp_clusters(TFlist, mzvalues):
    clusters = []
    curr_cluster = []
    for i, curr_point in enumerate(TFlist):
        if curr_point:
            curr_cluster.append(mzvalues[i])
        else:
            if curr_cluster:
                clusters.append(curr_cluster)
                curr_cluster = []
    return clusters


def find_nearest(array, value):
    """Index of the element in ``array`` closest to ``value``."""
    nearest_idx = np.where(abs(array - value) == abs(array - value).min())[0]
    return nearest_idx[0]


def get_unk_comp_clusters(ref_clusters, mzvalues):
    clusters = []
    mzvalues_array = np.array(mzvalues)
    mzvalues_array_float = np.asarray(mzvalues_array, dtype=np.float64, order='C')
    for i, curr_cluster in enumerate(ref_clusters):
        ref_ini_point = curr_cluster[0]
        indexnearest_unk = find_nearest(mzvalues_array_float, float(ref_ini_point))
        if ((indexnearest_unk + len(curr_cluster)) <= len(mzvalues)):
            clusters.append(mzvalues[range(indexnearest_unk, indexnearest_unk + len(curr_cluster))])
        curr_cluster = []
    return clusters


def spectrum_save(mz, intensity, name):
    originaltable = pd.DataFrame(np.array(mz))
    originaltable = originaltable.rename(columns={0: "mz"})
    originaltable["intensity"] = intensity
    originaltable.to_csv("{name}.csv".format(name=name))


def GKernal(x0, x, y, sig=1.):
    """Gaussian kernel smoother evaluated at ``x0``."""
    ax = x - x0
    Kernal = np.exp(-0.5 * np.square(ax) / np.square(sig))
    Weighted = Kernal * y
    return np.sum(Weighted) / np.sum(Kernal)


def gridding(mz, intensity, start=None, end=None, increment=0.125, epsilon=0.025):
    """Re-grid an irregular m/z axis onto a uniform grid with Nadaraya-Watson smoothing."""
    if start is None:
        start = math.floor(np.min(mz)) - 1
    if end is None:
        end = math.ceil(np.max(mz)) + 1
    target_grid = np.arange(start, end, increment)
    intensity_out = np.zeros(target_grid.shape[0])
    for i in range(target_grid.shape[0]):
        calculation = GKernal(target_grid[i], mz, intensity, sig=epsilon)
        if calculation == calculation:
            intensity_out[i] = calculation
        else:
            intensity_out[i] = 0
    return [target_grid, intensity_out]
