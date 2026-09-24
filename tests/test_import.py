"""Import-level smoke checks."""

import pgalign as pg


def test_version_is_string():
    assert isinstance(pg.__version__, str)
    assert pg.__version__.count(".") >= 1


def test_public_api_present():
    for name in [
        "AnnDataMALDI",
        "PeakCalling",
        "PeakCalling_single",
        "PeakCallingmv",
        "MALDI_SIM",
        "PGmzalign",
        "gridding",
        "GKernal",
    ]:
        assert hasattr(pg, name), f"pgalign is missing public name: {name}"


def test_paper_aligned_methods_exist():
    """Method names from the current PGAlign manuscript are exposed on the right classes."""
    assert hasattr(pg.PeakCalling, "peak_calling")
    assert hasattr(pg.PeakCalling, "peak_grouping")
    assert hasattr(pg.AnnDataMALDI, "peak_group_pairing")
    assert hasattr(pg.AnnDataMALDI, "fine_alignment_assessment")
