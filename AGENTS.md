# AGENTS.md

Guidance for coding agents working in this repository. Users of the package
should start from `README.md` or https://stagill.github.io/PGAlign/.

## What this is

PGAlign aligns an unknown mass spectrometry dataset onto a reference by peak
groups and returns the unknown spectra on the reference m/z grid. The package is
`pgalign` (on PyPI). The method is described in Deng, Zhang and Zhang,
"PGAlign: peak group alignment of mass spectrometry data", under review at the
*Journal of Proteome Research*; release v1.0.1 accompanies the manuscript.

## Set up and test

```bash
pip install -e ".[test,imzml]"
pytest -q
pip install -r docs/requirements.txt && mkdocs build --strict
```

CI (`.github/workflows/tests.yml`) runs the suite on Python 3.9 at the
dependency floors and on Python 3.14 with the latest releases. The floors are
stated in three places: `pyproject.toml`, the pins of the minimum leg in
`tests.yml`, and `docs/install.md`. Change all three together.
`tests/test_docs.py` runs the Python examples of `README.md` and `docs/`, so an
example that stops working fails the build.

## Layout

- `pgalign/pipeline.py`: `align()` and `AlignmentResult`, the one-call interface, and its input checks.
- `pgalign/PeakCalling.py`: Steps 1 and 2 (`PeakCallingmv`, `PeakCalling`, `PeakCalling_single`).
- `pgalign/AnnDataMALDI.py`: Steps 3 and 4 (`AnnDataMALDI`), the aligned dataset (`PGmzalign`) and the simulation (`MALDI_SIM`).
- `pgalign/readers.py`: `read_imzml`, `read_scils_csv`, `read_h5ad`, `write_h5ad`.
- `pgalign/util.py`: `gridding`, `save_clusterresults` and `spectrum_save`, plus helpers that serve the step classes and are not exported from `pgalign`.
- `pgalign/cli.py`: the `pgalign align` command.
- `pgalign/datasets.py` and `pgalign/data/`: the toy dataset; `scripts/make_toy_data.py` rebuilds it.
- `CodeInPaper/`: the notebooks behind the manuscript's analyses. They call the step classes directly.

## Conventions

- In memory, m/z values live in `var["m/z"]`. anndata 0.13 cannot write a
  column named `m/z`, so `pg.write_h5ad` stores it as `mz` and `pg.read_h5ad`
  restores it. Use these two rather than `write_h5ad` and `sc.read_h5ad`.
- The aligned dataset always lies on the reference grid and is cut to the
  shorter of the two m/z axes.
- `pg.align` names its parameters after the paper's symbols, while the step
  methods keep their historical names:

  | In the paper | `pg.align` | Step method and argument |
  |---|---|---|
  | &alpha;, peak calling quantile | `alpha` | `peak_calling(threshold=...)` |
  | peak grouping parameter | `percentile` | `peak_grouping(percentile=...)` |
  | &delta;, pairing criterion | `delta` | `peak_group_pairing(criteria=...)` |
  | &tau;, tolerance of the ambiguity diagnostic | `tau` | `flag_ambiguous_alignments(delta=...)` |

  `delta` therefore means &delta; in `pg.align` but &tau; in `flag_ambiguous_alignments`.
- `smooth=True`, the default, calls peaks with `PeakCallingmv` on the three-point
  moving average of the Methods, as in the atherosclerosis analysis;
  `smooth=False` uses `PeakCalling`, as in the simulation study and the canine
  sarcoma analysis. Every analysis in the manuscript ran Step 4 with
  `ignore=True`.
- The pre-2026 names (`callpeak`, `grouppeaks`, `greedy_match`, `fine_align`,
  and the keyword `threshould=`) remain as aliases that emit `DeprecationWarning`.

## Do not change the results silently

A change to `pgalign/` that alters a numerical output (called peaks, peak
groups, pairing, offsets, the aligned matrix, the ambiguity table) changes what
the manuscript reports. Before such a change, record these outputs on the toy
data for both peak callers; afterwards, show that they are identical, or explain
each difference in the commit message. Discrete outputs must match exactly.
Correlations may differ in the last digits across numpy versions (by up to
3e-7 between numpy 1.26 and 2.x).

## Releases

Publishing a GitHub release runs `.github/workflows/publish.yml`, which checks
that the tag equals `v` followed by `pgalign.__version__` and uploads the
package to PyPI. A version number can be uploaded only once. Bump
`__version__` in `pgalign/__init__.py` and `version` in `CITATION.cff` together.

## Style

Docstrings use the `:param:` style that mkdocstrings renders. Prose and
comments use no em-dashes.
