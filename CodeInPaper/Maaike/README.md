# Maaike — atherosclerosis regression macrophage metabolomics

Reproduces the macrophage metabolomics case study (Section 3.2 of the paper):
GALAXY aligns Week 5 MALDI spectra to Week 2, with subsequent joint spatial
segmentation and classification.

## Notebooks

| File                      | Purpose                                                                 |
|---------------------------|-------------------------------------------------------------------------|
| `Maaike_figure_a.ipynb`   | Main GALAXY alignment + Pearson correlation matrices + anchor points    |
| `MSIWarp.ipynb`           | Comparative analysis with MSIWarp on the same data                      |
| `Maaike_figure_b.Rmd`     | Spatial visualisation / clustering plots in R                           |
| `Maaike_figure_c_d.Rmd`   | Additional R-side figures                                               |

## Data

The MALDI data come from the atherosclerosis regression mouse model
(macrophage metabolomics, DHB positive and 9AA negative ionisation modes).
**The data are available from the data owners upon reasonable request** —
they are not redistributed in this repository.

Each notebook expects four CSV files per time point (Week 2 and Week 5)
named according to the original convention:

```
{DATA_DIR}/Maaike/{weekpoint} - {group} - All Spectra.csv
{DATA_DIR}/Maaike/{weekpoint} - {group} - Region Spots.csv
```

where `weekpoint` ∈ {`"2 wk regression"`, `"5 wk regression"`} and `group` ∈
{`"DHB pos"`, `"9AA neg"`}.

## How to run

1. Place the CSVs under `data/Maaike/` (or wherever you choose).
2. Open the notebook; edit the **CONFIGURE ME** cell so `DATA_DIR` points to
   the parent of the `Maaike/` folder.
3. Run the notebook end-to-end.
