# Tutorial

A minimal end-to-end walk-through. For a worked example on the atherosclerosis
data of the manuscript, see `Tutorial_PGAlign.ipynb` in the repository root.

## 1. Inputs

PGAlign operates on two `AnnData` objects: the *reference* spectrum and the
*unknown* spectrum to be aligned. Each `AnnData.X` is a (spectra &times; m/z bins)
intensity matrix; `var` carries m/z values as the `"m/z"` column; `obs` holds
per-pixel coordinates.

To read your own files, `pg.read_imzml` handles imzML and `pg.read_scils_csv`
SCiLS Lab CSV exports (see [Readers](api/readers.md)). You can also use the
bundled [toy dataset](toy-data.md) for a quick try-out:

```python
import pgalign as pg
unk, ref = pg.datasets.load_mouse_pancreas_toy()
```

## 2. Align

```python
res = pg.align(unk, ref)
print(res.summary())      # 7 groups paired, 7 shifted, 0 flagged at tau = 0.1
```

`pg.align` scales each spectrum to the median total intensity (Methods, "Data
normalization"), runs the four steps of the pipeline and the diagnostic for
ambiguous alignments, and leaves `unk` and `ref` unchanged. The defaults are
those of the Methods and of the atherosclerosis analysis in the manuscript:
peaks are called on the three-point moving average with &alpha; = 0.9
(`alpha`), the peak grouping parameter is 0.9 (`percentile`), the pairing
criterion is &delta; = 0 (`delta`), and the diagnostic uses &tau; = 0.1 (`tau`).
Pass `smooth=False` for the peak caller used in the simulation study and the
canine sarcoma analysis.

## 3. The results

```python
aligned = res.aligned       # the unknown spectra on the reference m/z grid
reference = res.reference   # the reference spectra on the same bins and scale
res.shifts                  # one row per paired peak group, with its offset
res.ambiguity               # shifted groups and the score of leaving them in place
```

`aligned` and `reference` share their m/z bins and their normalized scale, so
they can be stacked and passed to PCA, Harmony, clustering or classification.
They are cut to the shorter of the two m/z axes, so they can lack a few bins at
the top of the range; on the toy data they keep 599 of the reference's 601.
Save them with `pg.write_h5ad(aligned, "aligned.h5ad")`: anndata 0.13 and later cannot write a column named `m/z`, so `pg.write_h5ad` stores
it as `mz` and `pg.read_h5ad` restores it.

`res.shifts` gives, for each paired peak group, the reference m/z range it
occupies and its offset in grid points. `res.ambiguity` lists every shifted
group with the score of the selected alignment and of leaving the group where
it was; a group is flagged when the two lie within &tau; of each other, and such
groups deserve a look before downstream analysis.

## 4. Step by step

The notebooks under `CodeInPaper/` run the four steps one at a time. The
intermediate objects are also available from `pg.align` as `res.peaks` (Steps 1
and 2) and `res.steps` (Steps 3 and 4).

```python
import scanpy as sc
sc.pp.normalize_total(unk)
sc.pp.normalize_total(ref)

pc = pg.PeakCallingmv(unk, ref)      # pg.PeakCalling for smooth=False
pc.peak_calling(threshold=0.9)       # Step 1: alpha = 0.9 quantile
pc.peak_grouping(percentile=0.9)     # Step 2: peak grouping parameter

align = pg.AnnDataMALDI(unk, ref)
align.get_corr_peakgroup_refined(pc.jointcluster)   # Pearson similarity matrix
align.peak_group_pairing(criteria=0)                # Step 3: greedy pairing
align.fine_alignment_assessment(ignore=True)        # Step 4: rigid shift
align.summarize()

aligned = pg.PGmzalign(align).getAnnSim()
```

`pc.jointcluster` is the list of joint peak groups used as the alignment
scaffold. After `summarize()`:

- `align.unknownalign`: aligned m/z indices in the unknown spectrum
- `align.referenalign`: corresponding m/z indices in the reference
- `align.changerecord`: per-group integer-index shift applied

This sequence gives the same result as `pg.align(unk, ref)`. Unlike `pg.align`,
it normalizes `unk` and `ref` in place.

See [Case studies](case-studies.md) for full reproductions of the manuscript's
analyses.

## Notes on terminology

The current manuscript revision renamed the four pipeline steps. Pre-2026 names
(`callpeak`, `grouppeaks`, `greedy_match`, `fine_align`) are kept as
deprecated aliases and emit `DeprecationWarning`:

| New (paper)                   | Old (deprecated) |
|-------------------------------|------------------|
| `peak_calling`                | `callpeak`       |
| `peak_grouping`               | `grouppeaks`     |
| `peak_group_pairing`          | `greedy_match`   |
| `fine_alignment_assessment`   | `fine_align`     |
