"""Output step, moving-average peak caller and simulation on the bundled toy slice."""

import random

import numpy as np
import pytest

import pgalign as pg


@pytest.fixture(scope="module", params=["PeakCalling", "PeakCallingmv"])
def toy_alignment(request):
    unk, ref = pg.datasets.load_mouse_pancreas_toy()
    pc = getattr(pg, request.param)(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    align = pg.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    align.peak_group_pairing(criteria=0)
    align.fine_alignment_assessment(threshold=0.0, ignore=True)
    align.summarize()
    return unk, ref, align


def test_every_paired_group_recovers_the_planted_shift(toy_alignment):
    _, _, align = toy_alignment
    assert len(align.changerecord) > 0
    assert set(align.changerecord) == {2}


def test_output_step_puts_unknown_on_reference_grid(toy_alignment):
    unk, ref, align = toy_alignment
    aligned = pg.PGmzalign(align).getAnnSim()
    n = aligned.n_vars
    # the +2 shift pushes the last two bins past the end of the reference axis
    assert n == ref.n_vars - 2
    assert list(aligned.obs_names) == list(unk.obs_names)
    ref_mz = ref.var["m/z"].to_numpy(dtype=float)
    np.testing.assert_array_equal(aligned.var["m/z"].to_numpy(dtype=float), ref_mz[:n])
    mean_ref = np.asarray(ref.X).mean(axis=0)[:n]
    before = np.corrcoef(np.asarray(unk.X).mean(axis=0)[:n], mean_ref)[0, 1]
    after = np.corrcoef(np.asarray(aligned.X).mean(axis=0), mean_ref)[0, 1]
    assert before < 0.5
    assert after > 0.99


@pytest.mark.parametrize("shuffle, nregion", [(False, 1), (True, 1), (True, 2)])
def test_simulation_applies_the_recorded_shifts(toy_alignment, shuffle, nregion):
    _, ref, align = toy_alignment
    random.seed(0)
    np.random.seed(0)
    sim = pg.MALDI_SIM(align)
    data = sim.getAnnSim(shuffle=shuffle, sigma=0.1, nregion=nregion)
    shifts = sim.shiftdatadf
    if shuffle:
        assert sorted(shifts["currentshift_fromprev"]) == sorted(shifts["shift_fromprev"])
    else:
        assert list(shifts["currentshift_fromorg"]) == list(align.changerecord)
    # each inserted or deleted column lengthens or shortens the m/z axis by one bin
    assert data.n_vars == ref.n_vars + int(shifts["currentshift_fromprev"].sum())
    assert data.n_obs == ref.n_obs
    mz = data.var["m/z"].to_numpy(dtype=float)
    ref_mz = ref.var["m/z"].to_numpy(dtype=float)
    assert mz[0] == ref_mz[0]
    np.testing.assert_allclose(np.diff(mz), ref_mz[1] - ref_mz[0])
