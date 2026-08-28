"""Import-level smoke checks."""

import GalaxyPython as gx


def test_version_is_string():
    assert isinstance(gx.__version__, str)
    assert gx.__version__.count(".") >= 1


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
        assert hasattr(gx, name), f"GalaxyPython is missing public name: {name}"


def test_paper_aligned_methods_exist():
    """Method names from the current GALAXY manuscript are exposed on the right classes."""
    assert hasattr(gx.PeakCalling, "peak_calling")
    assert hasattr(gx.PeakCalling, "peak_grouping")
    assert hasattr(gx.AnnDataMALDI, "peak_group_pairing")
    assert hasattr(gx.AnnDataMALDI, "fine_alignment_assessment")
