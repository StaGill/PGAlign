"""Smoke test for PGAlign Steps 1 and 2 on synthetic data."""

import numpy as np

import pgalign as pg


def test_peak_calling_and_grouping_produce_clusters(tiny_anndata):
    unk, ref = tiny_anndata
    pc = pg.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    # At least one joint peak group must be detected on the four-peak fixture.
    assert len(pc.jointcluster) > 0
    # Each joint peak group is a list of m/z values.
    for group in pc.jointcluster:
        assert len(group) >= 1


def test_importpeak_keeps_the_two_spectra_apart(tiny_anndata):
    unk, ref = tiny_anndata
    pc = pg.PeakCalling(unk, ref)
    pc.importpeak(np.array([10, 20]), np.array([30, 40]))
    assert np.flatnonzero(pc.peakUnk).tolist() == [10, 20]
    assert np.flatnonzero(pc.peakRef).tolist() == [30, 40]


def test_imported_peaks_give_the_same_groups_as_step_1(tiny_anndata):
    unk, ref = tiny_anndata
    called = pg.PeakCalling(unk, ref)
    called.peak_calling(threshold=0.9)
    called.peak_grouping(percentile=0.9)
    imported = pg.PeakCalling(unk, ref)
    imported.importpeak(np.flatnonzero(called.peakUnk), np.flatnonzero(called.peakRef))
    imported.peak_grouping(percentile=0.9)
    assert imported.jointcluster == called.jointcluster


def test_single_spectrum_importpeak(tiny_anndata):
    _, ref = tiny_anndata
    single = pg.PeakCalling_single(ref)
    single.importpeak(None, np.array([30, 40]))
    assert np.flatnonzero(single.peakRef).tolist() == [30, 40]
