"""Import-level smoke checks."""

import pgalign as pg


def test_version_is_string():
    assert isinstance(pg.__version__, str)
    assert pg.__version__.count(".") >= 1


def test_public_api_present():
    for name in [
        "align",
        "AlignmentResult",
        "AnnDataMALDI",
        "PeakCalling",
        "PeakCalling_single",
        "PeakCallingmv",
        "MALDI_SIM",
        "PGmzalign",
        "gridding",
        "save_clusterresults",
        "spectrum_save",
        "read_imzml",
        "read_scils_csv",
        "read_h5ad",
        "write_h5ad",
    ]:
        assert hasattr(pg, name), f"pgalign is missing public name: {name}"


def test_removed_names_are_gone():
    """Single-cell leftovers were deleted in 2.0; internal helpers stay in pgalign.util only."""
    for name in ["prefilter_cells", "prefilter_genes", "prefilter_specialgenes", "GKernal",
                 "rotate", "group_range", "in_range_lookup", "comp_clusters", "find_nearest",
                 "get_unk_comp_clusters"]:
        assert not hasattr(pg, name), f"pgalign still exports {name}"
    assert callable(pg.util.in_range_lookup) and callable(pg.util.nadaraya_watson)


def test_paper_aligned_methods_exist():
    """Method names from the current PGAlign manuscript are exposed on the right classes."""
    assert hasattr(pg.PeakCalling, "peak_calling")
    assert hasattr(pg.PeakCalling, "peak_grouping")
    assert hasattr(pg.AnnDataMALDI, "peak_group_pairing")
    assert hasattr(pg.AnnDataMALDI, "fine_alignment_assessment")
