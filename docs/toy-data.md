# Toy data

A bundled 120-spectrum / 601-m/z-bin slice of the public mouse pancreas MALDI
dataset, useful for quick try-outs and for the test suite.

## Source

Sliced from the Zenodo deposit
[`10.5281/zenodo.3607915`](https://doi.org/10.5281/zenodo.3607915) (MALDI FTICR
MS imaging data of pancreatic mouse tissue; Prade et al., 2020; CC BY 4.0), after
gridding to 0.25 m/z in `CodeInPaper/MousePancreate.ipynb`. The bundled slice covers:

- Spatial window: 120 spectra drawn without replacement (seed 0) from the 21&times;21 pixels around `(x=1636, y=863)`.
- m/z range: `[100, 250]`, 601 bins at 0.25 m/z spacing.
- The "unknown" companion reads every spectrum two bins lower in m/z, so each
  peak sits two bins (0.5 m/z) higher than in the reference, and multiplies it
  by lognormal noise (&sigma; = 0.05, about &plusmn;5%). PGAlign should recover
  a per-group shift of +2.

The slice is stored as a compressed `.npz` (~230 KB) inside the package.
`scripts/make_toy_data.py` rebuilds it from the gridded source file.

## Usage

```python
import pgalign as pg

unk, ref = pg.datasets.load_mouse_pancreas_toy()

print(unk.shape, ref.shape)        # (120, 601) (120, 601)
print(ref.var["m/z"].head())       # 100.0, 100.25, ..., 250.0
```

## End-to-end on the toy data

```python
import pgalign as pg

unk, ref = pg.datasets.load_mouse_pancreas_toy()

pc = pg.PeakCalling(unk, ref)
pc.peak_calling(threshold=0.9)
pc.peak_grouping(percentile=0.9)

align = pg.AnnDataMALDI(unk, ref)
align.get_corr_peakgroup_refined(pc.jointcluster)
align.peak_group_pairing(criteria=0)
align.fine_alignment_assessment(threshold=0.0, ignore=True)
align.summarize()

import numpy as np
print("recovered mean shift:", float(np.mean(align.changerecord)))  # ~+2
```

## API

::: pgalign.datasets.load_mouse_pancreas_toy
