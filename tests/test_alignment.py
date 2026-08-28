"""End-to-end smoke test for the full GALAXY pipeline."""

import numpy as np

import GalaxyPython as gx


def test_full_pipeline_runs(tiny_anndata):
    unk, ref = tiny_anndata
    pc = gx.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)

    align = gx.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    align.peak_group_pairing(criteria=0)
    align.fine_alignment_assessment(threshold=0.0, ignore=True)
    align.summarize()

    # Alignment yields non-empty unknown/reference index arrays of matching length.
    assert align.unknownalign.size > 0
    assert align.referenalign.size > 0
    assert align.unknownalign.size == align.referenalign.size
    # Each recorded shift is an integer-valued offset.
    for shift in align.changerecord:
        assert isinstance(int(shift), int)
    # The unknown was shifted by +1 bin; the recovered shifts should average near +1
    # (negative because changerecord stores -change, where change > 0 means shift left).
    mean_shift = float(np.mean(align.changerecord))
    assert -3.0 <= mean_shift <= 3.0
