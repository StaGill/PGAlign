"""End-to-end alignment on the bundled mouse-pancreas toy slice."""

import numpy as np

import GalaxyPython as gx


def test_toy_dataset_round_trip():
    unk, ref = gx.datasets.load_mouse_pancreas_toy()
    assert unk.shape == ref.shape
    assert "m/z" in unk.var.columns

    pc = gx.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    assert len(pc.jointcluster) > 0

    align = gx.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    align.peak_group_pairing(criteria=0)
    align.fine_alignment_assessment(threshold=0.0, ignore=True)
    align.summarize()
    assert align.unknownalign.size > 0
    # The toy unknown is shifted by +2 bins; recovered shifts should average near +2.
    mean_shift = float(np.mean(align.changerecord))
    assert 0.5 <= mean_shift <= 3.5
