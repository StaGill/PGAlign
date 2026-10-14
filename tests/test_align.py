"""The one-call pipeline pgalign.align on the bundled toy slice."""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import pgalign as pg


@pytest.fixture(scope="module")
def toy():
    return pg.datasets.load_mouse_pancreas_toy()


@pytest.mark.parametrize("smooth, caller", [(True, "PeakCallingmv"), (False, "PeakCalling")])
def test_align_matches_the_step_by_step_pipeline(toy, smooth, caller):
    unk, ref = toy
    res = pg.align(unk, ref, smooth=smooth)
    unk, ref = unk.copy(), ref.copy()
    sc.pp.normalize_total(unk)
    sc.pp.normalize_total(ref)
    pc = getattr(pg, caller)(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    steps = pg.AnnDataMALDI(unk, ref)
    steps.get_corr_peakgroup_refined(pc.jointcluster)
    steps.peak_group_pairing(criteria=0)
    steps.fine_alignment_assessment(threshold=0.2, ignore=True)
    steps.summarize()
    aligned = pg.PGmzalign(steps).getAnnSim()
    np.testing.assert_array_equal(np.asarray(res.aligned.X), np.asarray(aligned.X))
    np.testing.assert_array_equal(res.aligned.var["m/z"].to_numpy(), aligned.var["m/z"].to_numpy())
    np.testing.assert_array_equal(res.steps.align_group, steps.align_group)
    assert list(res.shifts["offset"]) == list(steps.changerecord)
    pd.testing.assert_frame_equal(res.ambiguity, steps.flag_ambiguous_alignments(delta=0.1))


@pytest.mark.parametrize("smooth, n_groups", [(True, 7), (False, 15)])
def test_result_on_toy_data(toy, smooth, n_groups):
    unk, ref = toy
    res = pg.align(unk, ref, smooth=smooth)
    assert res.summary() == f"{n_groups} groups paired, {n_groups} shifted, 0 flagged at tau = 0.1"
    assert set(res.shifts["offset"]) == {2}
    assert (res.shifts["mz_start"] <= res.shifts["mz_end"]).all()
    assert res.aligned.shape == (unk.n_obs, ref.n_vars - 2)


def test_align_leaves_inputs_unchanged(toy):
    unk, ref = toy
    before = [np.asarray(unk.X).copy(), np.asarray(ref.X).copy()]
    pg.align(unk, ref)
    np.testing.assert_array_equal(np.asarray(unk.X), before[0])
    np.testing.assert_array_equal(np.asarray(ref.X), before[1])
