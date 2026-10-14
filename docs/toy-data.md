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
res = pg.align(unk, ref)

print(res.summary())                  # 7 groups paired, 7 shifted, 0 flagged at tau = 0.1
print(res.shifts["offset"].unique())  # [2]: every paired group recovers the planted shift

aligned = res.aligned                 # unknown spectra on the reference m/z grid
print(aligned.shape)                  # (120, 599): two bins lost at the top
```

With `smooth=False` (the peak caller of the simulation study and the canine
sarcoma analysis) 15 groups are paired, and each again recovers the shift of +2.

## API

::: pgalign.datasets.load_mouse_pancreas_toy
