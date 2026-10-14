# Tutorial

A minimal end-to-end walk-through. For a fully worked example with plots, see
`Tutorial_PGAlign.ipynb` in the repository root.

## 1. Inputs

PGAlign operates on two `AnnData` objects: the *reference* spectrum and the
*unknown* spectrum to be aligned. Each `AnnData.X` is a (spectra &times; m/z bins)
intensity matrix; `var` carries m/z values as the `"m/z"` column; `obs` holds
per-pixel coordinates.

You can use the bundled [toy dataset](toy-data.md) for a quick try-out:

```python
import pgalign as pg
unk, ref = pg.datasets.load_mouse_pancreas_toy()
```

## 2. Per-spectrum normalisation

```python
import scanpy as sc
sc.pp.normalize_per_cell(unk)
sc.pp.normalize_per_cell(ref)
```

## 3. Peak calling and peak grouping (Steps 1 & 2)

```python
pc = pg.PeakCalling(unk, ref)
pc.peak_calling(threshold=0.9)        # Step 1: alpha = 0.9 quantile
pc.peak_grouping(percentile=0.9)      # Step 2: top 10% inter-peak distance cut
```

`pc.jointcluster` is the list of joint peak groups used as the alignment scaffold.

## 4. Peak group pairing and fine alignment (Steps 3 & 4)

```python
align = pg.AnnDataMALDI(unk, ref)
align.get_corr_peakgroup_refined(pc.jointcluster)   # Pearson similarity matrix
align.peak_group_pairing(criteria=0)                # Step 3: greedy pairing
align.fine_alignment_assessment(threshold=0.2, ignore=True)  # Step 4: rigid shift
align.summarize()
```

After `summarize()`:

- `align.unknownalign`: aligned m/z indices in the unknown spectrum
- `align.referenalign`: corresponding m/z indices in the reference
- `align.changerecord`: per-group integer-index shift applied

## 5. The aligned dataset

```python
aligned = pg.PGmzalign(align).getAnnSim()
```

`aligned` holds the unknown spectra on the reference m/z grid. Its
`var["m/z"]` are the first `aligned.n_vars` m/z values of the reference, so it
can be stacked with `ref[:, :aligned.n_vars]` and passed to PCA, Harmony,
clustering or classification. Despite the method's name, nothing is simulated.
The result is cut to the shorter of the two m/z axes, so it can lack a few bins
at the top of the range; on the toy data it keeps 599 of the reference's 601.

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
