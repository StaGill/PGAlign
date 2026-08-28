# Canine sarcomas — alignment and classification

Reproduces the canine sarcoma case study (Section 3.3 of the paper). GALAXY
aligns cancer-tissue spectra (samples Ds1–Ds4, Ds18–Ds20, Ds24) to a normal
tissue reference (Ds26); classification performance is then evaluated before
and after alignment, alongside joint clustering across cancer samples.

## Notebooks

| File                                 | Purpose                                                                       |
|--------------------------------------|-------------------------------------------------------------------------------|
| `Sarcomas_Analysis.ipynb`            | PCA + LDA classification (leave-one-out CV) on `.h5ad` files for each sample  |
| `Sarcomas_figure_a.ipynb`            | Figure-A generation: alignment heatmaps and shift distributions               |
| `Sarcomas_figure_b_c_d.Rmd`          | Figures B–D in R                                                              |
| `Clustering (sarcomas)/Harmony_ds2_24.Rmd` | Harmony-based joint clustering of two samples after GALAXY alignment   |

## Data

Source: ProteomeXchange Consortium (PRIDE) accession `PXD010990`
(<https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD010990>).

The notebooks expect pre-processed `AnnData` files at:

```
{DATA_DIR}/CanineSarcomas/CanineSarcomas_Ds{ID}.h5ad        # unaligned (default gridding)
{DATA_DIR}/CanineSarcomas/CanineSarcomas_Ds{ID}_05.h5ad     # gridded at 0.5 m/z
```

where `ID` ∈ {`1`, `2`, `3`, `4`, `18`, `19`, `20`, `24`, `26`}.

Preprocessing from raw PRIDE files to `.h5ad` follows the standard pipeline:
mass-range truncation, gridding (Nadaraya–Watson, `increment=0.5`), TIC
normalisation, baseline filtering (> 2.2 × 85th percentile), spatial-frequency
filtering (< 0.5% occurrence removed). See Methods §Data Preprocessing in the
paper.

## How to run

1. Download `PXD010990` from PRIDE and preprocess each sample to `.h5ad`
   following the steps above (or use your own pipeline producing the same
   layout).
2. Place the `.h5ad` files under `data/CanineSarcomas/`.
3. Open the notebook; edit the **CONFIGURE ME** cell so `DATA_DIR` points to
   the parent of `CanineSarcomas/`.
4. Run end-to-end.
