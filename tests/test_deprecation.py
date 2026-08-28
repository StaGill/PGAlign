"""Deprecated aliases still work and emit DeprecationWarning."""

import warnings

import pytest

import GalaxyPython as gx


def test_callpeak_alias_warns(tiny_anndata):
    unk, ref = tiny_anndata
    pc = gx.PeakCalling(unk, ref)
    with pytest.warns(DeprecationWarning, match="peak_calling"):
        pc.callpeak(threshold=0.9)


def test_grouppeaks_alias_warns(tiny_anndata):
    unk, ref = tiny_anndata
    pc = gx.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    with pytest.warns(DeprecationWarning, match="peak_grouping"):
        pc.grouppeaks(percentile=0.9)


def test_greedy_match_alias_warns(tiny_anndata):
    unk, ref = tiny_anndata
    pc = gx.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    align = gx.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    with pytest.warns(DeprecationWarning, match="peak_group_pairing"):
        align.greedy_match(criteria=0)


def test_fine_align_alias_warns(tiny_anndata):
    unk, ref = tiny_anndata
    pc = gx.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    align = gx.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    align.peak_group_pairing(criteria=0)
    with pytest.warns(DeprecationWarning, match="fine_alignment_assessment"):
        align.fine_align(threshould=0.0, ignore=True)


def test_threshould_kwarg_alias_warns(tiny_anndata):
    """The fine_alignment_assessment(threshould=...) typo alias still works."""
    unk, ref = tiny_anndata
    pc = gx.PeakCalling(unk, ref)
    pc.peak_calling(threshold=0.9)
    pc.peak_grouping(percentile=0.9)
    align = gx.AnnDataMALDI(unk, ref)
    align.get_corr_peakgroup_refined(pc.jointcluster)
    align.peak_group_pairing(criteria=0)
    with pytest.warns(DeprecationWarning, match="threshold"):
        align.fine_alignment_assessment(threshould=0.0, ignore=True)
