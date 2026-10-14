"""One-call interface to the PGAlign pipeline."""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from .AnnDataMALDI import AnnDataMALDI, PGmzalign
from .PeakCalling import PeakCalling, PeakCallingmv


def _check_inputs(unk, ref):
    """Raise a ValueError that names the problem and the fix for input pg.align cannot use."""
    ranges = []
    for name, adata in (("unk", unk), ("ref", ref)):
        if "m/z" not in adata.var.columns:
            hint = ("; a file saved with pg.write_h5ad stores them as 'mz', so read it with "
                    "pg.read_h5ad") if "mz" in adata.var.columns else ""
            raise ValueError(
                f"{name}.var has no 'm/z' column (columns: {list(adata.var.columns)}); "
                f"store the m/z values as {name}.var['m/z']{hint}"
            )
        mz = adata.var["m/z"].to_numpy(dtype=float)
        if not np.isfinite(mz).all():
            raise ValueError(f"{name}.var['m/z'] contains missing or infinite values")
        if np.any(np.diff(mz) <= 0):
            raise ValueError(
                f"{name}.var['m/z'] must be strictly increasing; sort the bins with "
                f"{name} = {name}[:, np.argsort({name}.var['m/z'].to_numpy())]"
            )
        values = adata.X.data if sp.issparse(adata.X) else np.asarray(adata.X)
        n_bad = int(np.count_nonzero(~np.isfinite(values)))
        if n_bad:
            raise ValueError(
                f"{name}.X contains {n_bad} missing or infinite intensities; "
                "replace them, for example with 0"
            )
        mean = np.asarray(adata.X.mean(axis=0)).ravel()
        if np.ptp(mean) == 0:
            raise ValueError(
                f"{name} carries no signal: its mean spectrum is constant at {mean[0]:g}, "
                "so no peaks can be called"
            )
        ranges.append((mz[0], mz[-1]))
    (unk_lo, unk_hi), (ref_lo, ref_hi) = ranges
    if max(unk_lo, ref_lo) >= min(unk_hi, ref_hi):
        raise ValueError(
            f"the m/z ranges do not overlap: unk covers {unk_lo:g} to {unk_hi:g}, "
            f"ref covers {ref_lo:g} to {ref_hi:g}"
        )


class AlignmentResult(object):
    """The output of :func:`align`.

    :ivar aligned: the unknown spectra on the reference m/z grid, as an
        ``anndata.AnnData``. ``var["m/z"]`` holds the reference m/z values and
        ``obs`` the unknown spectra. It is cut to the shorter of the two m/z axes,
        so it can lack a few bins at the top of the range.
    :ivar reference: the reference spectra on the same bins and the same scale as
        ``aligned``, so the two can be analysed jointly. Both carry the
        reference's bin labels as ``var_names``.
    :ivar shifts: a DataFrame with one row per paired peak group: ``unknown_group``,
        ``reference_group``, ``mz_start`` and ``mz_end`` (reference m/z range) and
        ``offset``, the shift in grid points applied in Step 4.
    :ivar ambiguity: the table of ``AnnDataMALDI.flag_ambiguous_alignments``, one
        row per shifted group.
    :ivar peaks: the peak calling object of Steps 1 and 2.
    :ivar steps: the ``AnnDataMALDI`` object of Steps 3 and 4.
    :ivar tau: the tolerance the ambiguity table was computed with.
    """

    def __init__(self, aligned, reference, shifts, ambiguity, peaks, steps, tau):
        self.aligned = aligned
        self.reference = reference
        self.shifts = shifts
        self.ambiguity = ambiguity
        self.peaks = peaks
        self.steps = steps
        self.tau = tau

    def summary(self):
        """Return a one-line count of paired, shifted and flagged peak groups."""
        n_flagged = int(self.ambiguity["flagged"].sum())
        return (f"{len(self.shifts)} groups paired, {len(self.ambiguity)} shifted, "
                f"{n_flagged} flagged at tau = {self.tau}")

    def __repr__(self):
        return f"AlignmentResult: {self.summary()}"


def align(unk, ref, smooth=True, alpha=0.9, percentile=0.9, delta=0, tau=0.1, normalize=True):
    """Align an unknown dataset onto a reference with the full PGAlign pipeline.

    Runs the four steps of the Methods in order, then the diagnostic for ambiguous
    alignments, and returns the aligned dataset. The defaults are those of the
    Methods and of the atherosclerosis analysis in the manuscript. ``smooth=True``
    calls peaks on the three-point moving average described in the Methods
    ("Peak calling"), with ``PeakCallingmv``; ``smooth=False`` uses
    ``PeakCalling``, as in the simulation study and the canine sarcoma analysis.
    Step 4 applies the best offset to every paired group, as every analysis in the
    manuscript did (``fine_alignment_assessment(ignore=True)``). The inputs are not
    modified; a sparse intensity matrix is accepted and densified on the copy.
    Input that cannot be aligned (no ``var["m/z"]``, m/z not increasing, missing
    intensities, m/z ranges that do not overlap) raises a ``ValueError`` saying
    what to fix.

    :param unk: the unknown dataset, an ``anndata.AnnData`` with the m/z values in
        ``var["m/z"]``.
    :param ref: the reference dataset, in the same form.
    :param smooth: call peaks on the moving average (``PeakCallingmv``) or on the
        raw mean spectrum (``PeakCalling``). Defaults to True.
    :param alpha: the peak calling quantile (alpha in the paper). Defaults to 0.9.
    :param percentile: the peak grouping parameter. Defaults to 0.9.
    :param delta: the pairing criterion (delta in the paper). Defaults to 0.
    :param tau: the tolerance of the diagnostic for ambiguous alignments (tau in the
        paper). Defaults to 0.1.
    :param normalize: scale each spectrum to the median total intensity before
        aligning (Methods, "Data normalization"), on copies of the inputs.
        ``aligned`` and ``reference`` are returned on this normalized scale.
        Defaults to True.
    :return: an :class:`AlignmentResult`.
    """
    _check_inputs(unk, ref)
    unk = unk.copy()
    ref = ref.copy()
    for adata in (unk, ref):
        if sp.issparse(adata.X):
            adata.X = adata.X.toarray()
    if normalize:
        sc.pp.normalize_total(unk)
        sc.pp.normalize_total(ref)
    peaks = (PeakCallingmv if smooth else PeakCalling)(unk, ref)
    peaks.peak_calling(threshold=alpha)
    peaks.peak_grouping(percentile=percentile)
    steps = AnnDataMALDI(unk, ref)
    steps.get_corr_peakgroup_refined(peaks.jointcluster)
    steps.peak_group_pairing(criteria=delta)
    steps.fine_alignment_assessment(ignore=True)
    steps.summarize()
    ambiguity = steps.flag_ambiguous_alignments(delta=tau)
    aligned = PGmzalign(steps).getAnnSim()
    reference = ref[:, :aligned.n_vars].copy()
    aligned.var_names = reference.var_names
    mz_ref = ref.var["m/z"].to_numpy(dtype=float)
    shifts = pd.DataFrame({
        "unknown_group": steps.align_group[:, 0].astype(int),
        "reference_group": steps.align_group[:, 1].astype(int),
        "mz_start": [mz_ref[group[0]] for group in steps.aligned_mz_clusters_ref],
        "mz_end": [mz_ref[group[-1]] for group in steps.aligned_mz_clusters_ref],
        "offset": np.asarray(steps.changerecord, dtype=int),
    })
    return AlignmentResult(aligned, reference, shifts, ambiguity, peaks, steps, tau)
