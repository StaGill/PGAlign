"""Smoke test for GALAXY Steps 1 and 2 on synthetic data."""

import GalaxyPython as gx


def test_peak_calling_and_grouping_produce_clusters(tiny_anndata):
    unk, ref = tiny_anndata
    pc = gx.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    # At least one joint peak group must be detected on the four-peak fixture.
    assert len(pc.jointcluster) > 0
    # Each joint peak group is a list of m/z values.
    for group in pc.jointcluster:
        assert len(group) >= 1
