"""Diagnostic for ambiguous alignments (AnnDataMALDI.flag_ambiguous_alignments)."""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import pgalign as pg


def _controlled_align(pearson, blocks, group_len=5, gap=20):
    """AnnDataMALDI whose similarity scores are set by hand.

    ``pearson`` is the group-level matrix (rows unknown, columns reference).
    ``blocks[(u, r)]`` maps a diagonal ``k`` of the 9 x 9 sliding block (numpy
    orientation, rows are unknown slides) to the correlation placed on it.
    """
    n_groups = len(pearson)
    mz = np.arange(100.0, 100.0 + gap * (n_groups + 1), 1.0)
    var = pd.DataFrame({"m/z": mz}, index=[str(v) for v in mz])
    adata = sc.AnnData(X=np.ones((2, len(mz))), var=var)
    align = pg.AnnDataMALDI(adata, adata.copy())
    align.ref_clusters = [list(var.index[10 + gap * g:10 + gap * g + group_len]) for g in range(n_groups)]
    align.unk_clusters = align.ref_clusters
    align.midpoints_unk = [c[len(c) // 2] for c in align.ref_clusters]
    align.nclusters = n_groups
    align.PearsonMatrix = np.asarray(pearson, dtype=float)
    full = np.zeros((9 * n_groups, 9 * n_groups))
    for (u, r), diagonals in blocks.items():
        block = np.zeros((9, 9))
        for a in range(9):
            for b in range(9):
                block[a, b] = diagonals.get(b - a, 0.0)
        full[u * 9:u * 9 + 9, r * 9:r * 9 + 9] = block
    align.PearsonMatrixFull = full
    return align


def test_offset_flagged_when_staying_is_almost_as_good():
    align = _controlled_align([[0.9]], {(0, 0): {-1: 0.80, 0: 0.75}})
    align.peak_group_pairing()
    align.fine_alignment_assessment(threshold=0.2, ignore=True)
    tab = align.flag_ambiguous_alignments(delta=0.1)
    assert align.ambiguity_table is tab
    assert len(tab) == 1
    row = tab.iloc[0]
    assert row["rule"] == "offset"
    assert row["offset"] == 1 == align.changerecord[0]
    assert row["score_selected"] == pytest.approx(0.80)
    assert row["score_in_place"] == pytest.approx(0.75)
    assert row["margin"] == pytest.approx(0.05)
    assert bool(row["flagged"])


def test_offset_not_flagged_when_shift_is_clearly_better():
    align = _controlled_align([[0.9]], {(0, 0): {2: 0.9, 0: 0.3}})
    align.peak_group_pairing()
    align.fine_alignment_assessment(threshold=0.2, ignore=True)
    tab = align.flag_ambiguous_alignments(delta=0.1)
    assert len(tab) == 1
    assert tab["offset"].iloc[0] == -2
    assert tab["margin"].iloc[0] == pytest.approx(0.6)
    assert not tab["flagged"].iloc[0]
    assert align.flag_ambiguous_alignments(delta=0.7)["flagged"].iloc[0]


@pytest.mark.parametrize("in_place, flagged", [(0.75, True), (0.3, False), (-1.0, False)])
def test_group_rule_compares_with_in_place_similarity(in_place, flagged):
    pearson = [[0.9, 0.0, 0.0],
               [0.0, in_place, 0.8],
               [0.0, 0.0, 0.0]]
    align = _controlled_align(pearson, {(0, 0): {0: 0.9}, (1, 2): {0: 0.5}})
    align.peak_group_pairing()
    assert align.align_group.tolist() == [[0, 0], [1, 2]]
    align.fine_alignment_assessment(threshold=0.2, ignore=True)
    tab = align.flag_ambiguous_alignments(delta=0.1)
    # group 0 stays in place with a zero offset and is not reported
    assert len(tab) == 1
    row = tab.iloc[0]
    assert (row["unknown_group"], row["reference_group"], row["rule"]) == (1, 2, "group")
    assert row["score_selected"] == pytest.approx(0.8)
    if in_place == -1.0:
        assert np.isnan(row["score_in_place"])
    else:
        assert row["score_in_place"] == pytest.approx(in_place)
    assert bool(row["flagged"]) is flagged


def test_groups_skipped_by_threshold_keep_bookkeeping():
    pearson = np.diag([0.9, 0.9, 0.9])
    blocks = {
        (0, 0): {1: 0.1, 0: 0.05},     # best diagonal below threshold: not accepted
        (1, 1): {-2: 0.9, 0: 0.85},    # accepted, shifted, nearly tied with staying
        (2, 2): {0: 0.9},              # accepted, left in place
    }
    align = _controlled_align(pearson, blocks)
    align.peak_group_pairing()
    align.fine_alignment_assessment(threshold=0.2)
    assert align.changerecord == [2, 0]
    tab = align.flag_ambiguous_alignments(delta=0.1)
    assert len(tab) == 1
    row = tab.iloc[0]
    assert row["aligned_group"] == 0
    assert row["unknown_group"] == 1 and row["reference_group"] == 1
    assert row["offset"] == 2
    assert row["mz_start"] == float(align.ref_clusters[1][0])
    assert row["mz_end"] == float(align.ref_clusters[1][-1])
    assert bool(row["flagged"])


def _check_table(align, delta):
    tab = align.flag_ambiguous_alignments(delta=delta)
    change = np.asarray(align.changerecord)
    ag = align.align_group
    accepted_refs = [c[0] for c in align.aligned_mz_clusters_ref]
    n_shifted = 0
    for i in range(ag.shape[0]):
        head = list(align.mz_valueRef.index).index(align.ref_clusters[ag[i, 1]][0])
        if head in accepted_refs:
            j = accepted_refs.index(head)
            n_shifted += int(ag[i, 0] != ag[i, 1] or change[j] != 0)
    assert len(tab) == n_shifted
    assert (tab["offset"].to_numpy() == change[tab["aligned_group"].to_numpy()]).all()
    known = tab["margin"].notna()
    assert (tab.loc[known, "flagged"] == (tab.loc[known, "margin"] <= delta)).all()
    assert not tab.loc[~known, "flagged"].any()
    offset_rows = tab[tab["rule"] == "offset"]
    assert (offset_rows["score_selected"] >= offset_rows["score_in_place"]).all()
    return tab


def test_pipeline_on_synthetic_spectra(tiny_anndata):
    unk, ref = tiny_anndata
    pc = pg.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    align = pg.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    align.peak_group_pairing(criteria=0)
    align.fine_alignment_assessment(threshold=0.2, ignore=True)
    counts = [_check_table(align, d)["flagged"].sum() for d in (0.05, 0.1, 0.2)]
    assert counts == sorted(counts)


def test_toy_dataset_sanity():
    try:
        unk, ref = pg.datasets.load_mouse_pancreas_toy()
    except FileNotFoundError:
        pytest.skip("bundled toy dataset is not available in this checkout")
    pc = pg.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    align = pg.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    align.peak_group_pairing(criteria=0)
    align.fine_alignment_assessment(threshold=0.0, ignore=True)
    tab = _check_table(align, 0.1)
    # the toy unknown carries a rigid +2-bin shift, so shifted groups must be found
    assert len(tab) > 0
