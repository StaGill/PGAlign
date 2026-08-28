# Toy data

A bundled 120-spectrum / 601-m/z-bin slice of the public mouse pancreas MALDI
dataset, useful for quick try-outs and for the test suite.

## Source

Sliced from the Zenodo deposit
[`10.5281/zenodo.3607915`](https://doi.org/10.5281/zenodo.3607915) (Mouse
islets of Langerhans MALDI, FT-ICR negative mode). The bundled slice covers:

- Spatial window: 21&times;21 pixels around `(x=1636, y=863)`, subsampled to 120 spectra (seed = 0).
- m/z range: `[100, 250]` &mdash; 601 bins at 0.25 m/z spacing.
- The "unknown" companion was created by shifting the reference by +2 m/z bins
  and applying mild lognormal noise (&sigma; = 0.05), so GALAXY should recover
  a per-group shift of about +2.

The slice is stored as a compressed `.npz` (~230 KB) inside the package.

## Usage

```python
import GalaxyPython as gx

unk, ref = gx.datasets.load_mouse_pancreas_toy()

print(unk.shape, ref.shape)        # (120, 601) (120, 601)
print(ref.var["m/z"].head())       # 100.0, 100.25, ..., 250.0
```

## End-to-end on the toy data

```python
import GalaxyPython as gx

unk, ref = gx.datasets.load_mouse_pancreas_toy()

pc = gx.PeakCalling(unk, ref)
pc.peak_calling(threshold=0.9)
pc.peak_grouping(percentile=0.9)

align = gx.AnnDataMALDI(unk, ref)
align.get_corr_peakgroup_refined(pc.jointcluster)
align.peak_group_pairing(criteria=0)
align.fine_alignment_assessment(threshold=0.0, ignore=True)
align.summarize()

import numpy as np
print("recovered mean shift:", float(np.mean(align.changerecord)))  # ~+2
```

## API

::: GalaxyPython.datasets.load_mouse_pancreas_toy
