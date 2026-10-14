"""PGAlign: peak group alignment of mass spectrometry data.

Implementation accompanying the Deng, Zhang & Zhang (2026) manuscript (under review).
The four pipeline steps in the paper map to these public objects:

  Step 1  Peak Calling                  -> PeakCalling.peak_calling()
  Step 2  Peak Grouping                 -> PeakCalling.peak_grouping()
  Step 3  Peak Group Pairing            -> AnnDataMALDI.peak_group_pairing()
  Step 4  Fine Alignment Assessment     -> AnnDataMALDI.fine_alignment_assessment()

align() runs all four with the manuscript's defaults and returns the aligned dataset.
"""

__version__ = "1.1.0"

from .util import (
    prefilter_cells,
    prefilter_genes,
    prefilter_specialgenes,
    save_clusterresults,
    rotate,
    group_range,
    in_range_lookup,
    comp_clusters,
    find_nearest,
    get_unk_comp_clusters,
    spectrum_save,
    GKernal,
    gridding,
)
from .AnnDataMALDI import (
    AnnDataMALDI,
    MALDI_SIM,
    PGmzalign,
)
from .PeakCalling import (
    PeakCalling,
    PeakCalling_single,
    PeakCallingmv,
)
from .pipeline import (
    AlignmentResult,
    align,
)
from .readers import (
    read_h5ad,
    read_imzml,
    read_scils_csv,
    write_h5ad,
)
from . import datasets

__all__ = [
    "__version__",
    # One-call pipeline
    "align",
    "AlignmentResult",
    # Readers and .h5ad I/O
    "read_imzml",
    "read_scils_csv",
    "read_h5ad",
    "write_h5ad",
    # AnnDataMALDI core
    "AnnDataMALDI",
    "MALDI_SIM",
    "PGmzalign",
    # Peak calling
    "PeakCalling",
    "PeakCalling_single",
    "PeakCallingmv",
    # Toy data
    "datasets",
    # Utilities
    "prefilter_cells",
    "prefilter_genes",
    "prefilter_specialgenes",
    "save_clusterresults",
    "rotate",
    "group_range",
    "in_range_lookup",
    "comp_clusters",
    "find_nearest",
    "get_unk_comp_clusters",
    "spectrum_save",
    "GKernal",
    "gridding",
]
