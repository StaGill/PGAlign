# GALAXY

**GALAXY (Group Alignment of Mass Spectrometry data)** is a peak-group-based
algorithm for aligning mass spectrometry (MS) spectra onto a common m/z grid.
It is designed for imaging and spatial metabolomics data such as MALDI-MS,
where spectra from different runs or tissues often exhibit small m/z shifts
that prevent direct comparison.

Anji Deng, Yuyang Zhang, and Qihuang Zhang\*
McGill University — corresponding author: <qihuang.zhang@mcgill.ca>

**Documentation:** <https://stagill.github.io/GALAXY/>

![Overview of the GALAXY workflow](images/galaxy_workflow.png)

By aligning an "unknown" spectrum (or dataset) to a reference while *forcing*
matched spectra to share the same m/z values, GALAXY enables downstream
analyses such as:

- joint spatial segmentation across tissues or time points
- classification using combined datasets
- other multi-sample analyses that require a common m/z grid

This repository contains the reference Python implementation accompanying the
GALAXY manuscript (Deng, Zhang & Zhang, 2026; under peer review).

---

## Repository layout

- `GalaxyPython/` — Core Python implementation of GALAXY (alignment and peak-group functions).
- `CodeInPaper/` — Scripts that reproduce the figures and results in the manuscript. See `CodeInPaper/README.md`.
- `Tutorial_GALAXY.ipynb` — A Jupyter notebook that walks through aligning two MALDI datasets (e.g., Week 2 and Week 5 macrophage samples) and preparing them for joint segmentation.
- `tests/` — Smoke tests (`pytest`).
- `CITATION.cff` — Citation metadata.
- `LICENSE` — MIT license.

---

## Installation

Clone the repository and install with `pip`:

```bash
git clone https://github.com/StaGill/GALAXY.git
cd GALAXY
pip install -e .
```

The package's importable name is `GalaxyPython`:

```python
import GalaxyPython as gx
print(gx.__version__)
```

---

## API at a glance

The four-step pipeline in the paper maps to the following public methods:

| Step | Manuscript section          | Object · method                                   |
|-----:|-----------------------------|---------------------------------------------------|
| 1    | Peak Calling                | `PeakCalling.peak_calling(threshold=0.9)`         |
| 2    | Peak Grouping               | `PeakCalling.peak_grouping(percentile=0.9)`       |
| 3    | Peak Group Pairing          | `AnnDataMALDI.peak_group_pairing(criteria=0)`     |
| 4    | Fine Alignment Assessment   | `AnnDataMALDI.fine_alignment_assessment(threshold=0.2)` |

Similarity is measured with Pearson's correlation; Step 3 uses a sliding
window of +/- 4 m/z units (see Methods of the paper).

The pre-2026 method names (`callpeak`, `grouppeaks`, `greedy_match`,
`fine_align`) are kept as deprecated aliases and emit `DeprecationWarning`.

---

## Toy data

A 120-spectrum / 601-m/z-bin slice of the public mouse pancreas dataset ships
with the package for tutorials and tests:

```python
import GalaxyPython as gx
unk, ref = gx.datasets.load_mouse_pancreas_toy()
```

The "unknown" spectrum has a planted +2-bin rigid shift relative to the
reference, so GALAXY should recover a shift of about +2 m/z bins.

---

## Quick start

A minimal example is given in `Tutorial_GALAXY.ipynb`. The basic workflow is:

1. Load two MALDI datasets (e.g., two time points or two tissues) with optional spatial coordinates.
2. Wrap them as `AnnData` objects (via `scanpy` / `anndata`).
3. Normalize intensities per spectrum.
4. Run the four GALAXY steps.
5. Extract the aligned spectra for joint downstream analyses (e.g., PCA + Harmony + clustering).

```python
import scanpy as sc
import GalaxyPython as gx

# X = spectra, var = m/z table, obs = coordinates (see Tutorial_GALAXY.ipynb)
MALDIdataAnn1 = ...   # unknown spectrum
MALDIdataAnn2 = ...   # reference spectrum

# Per-spectrum normalisation
sc.pp.normalize_per_cell(MALDIdataAnn1)
sc.pp.normalize_per_cell(MALDIdataAnn2)

# Step 1 + Step 2: peak calling and peak grouping
PeakGroup = gx.PeakCalling(MALDIdataAnn1, MALDIdataAnn2)
PeakGroup.peak_calling(threshold=0.9)
PeakGroup.peak_grouping(percentile=0.9)

# Step 3: peak group pairing
ExactAlign = gx.AnnDataMALDI(MALDIdataAnn1, MALDIdataAnn2)
ExactAlign.get_corr_peakgroup_refined(PeakGroup.jointcluster)
ExactAlign.peak_group_pairing(criteria=0)

# Step 4: fine alignment assessment
ExactAlign.fine_alignment_assessment(threshold=0.2, ignore=True)
ExactAlign.summarize()
```

Please refer to `Tutorial_GALAXY.ipynb` for a fully worked example with real
data and plots.

---

## Reproducing results from the manuscript

The notebooks under `CodeInPaper/` reproduce the main analyses. Each
subfolder has its own `README.md` with dataset access instructions:

- **Simulation study (mouse pancreas MALDI)** — alignment error under
  controlled peak shifts and noise. See `CodeInPaper/MousePancreate.ipynb`.
- **Macrophage metabolomics (atherosclerosis regression)** — aligns Week 5
  MALDI spectra to Week 2 and performs joint spatial segmentation. See
  `CodeInPaper/Maaike/`.
- **Canine sarcoma classification** — aligns cancer tissue spectra to normal
  tissue and assesses classification performance before/after alignment. See
  `CodeInPaper/Sarcomas/`.

**Data sources:**

- Mouse pancreas MALDI: <https://doi.org/10.5281/zenodo.3607915>
- Atherosclerosis regression MALDI: available from the data owners upon reasonable request.
- Canine carcinoma MALDI: ProteomeXchange PRIDE accession `PXD010990` (<https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD010990>).

---

## Citation

The GALAXY manuscript is currently under peer review. If you use GALAXY in
your research, please cite this repository (see `CITATION.cff`, which GitHub
renders as a "Cite this repository" button) and contact the corresponding
author for the latest preprint:

> Deng A, Zhang Y, Zhang Q. *GALAXY: Group Alignment of Mass Spectrometry data
> for imaging and spatial metabolomics*. Manuscript under review, 2026.

---

## License

MIT — see `LICENSE`.
