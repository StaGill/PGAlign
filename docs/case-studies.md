# Case studies

The notebooks under [`CodeInPaper/`](https://github.com/StaGill/PGAlign/tree/main/CodeInPaper)
reproduce the analyses in the PGAlign manuscript.

| Subfolder / file                  | Case study                              | Results subsection             | Dataset access                                                                 |
|-----------------------------------|-----------------------------------------|--------------------------------|--------------------------------------------------------------------------------|
| `MousePancreate.ipynb`            | Simulation on mouse pancreas MALDI       | Simulation studies with mass spectrometry imaging data of mouse pancreatic tissue | Zenodo: <https://doi.org/10.5281/zenodo.3607915>                               |
| `Maaike/`                         | Atherosclerosis regression macrophages   | Application of PGAlign to MALDI data from an atherosclerosis regression model | Mean spectra on Zenodo: <https://zenodo.org/records/20533948>; per-spot spectra from the data owners upon reasonable request |
| `Sarcomas/`                       | Canine sarcoma classification             | Exploration of the ability to preserve the biological difference in dog sarcoma tumor data | ProteomeXchange PRIDE accession `PXD010990`                                    |

Each subfolder has its own `README.md` with dataset access details and a
**CONFIGURE ME** cell at the top of each notebook that sets `DATA_DIR`,
`OUTPUT_DIR`, and `REPO_ROOT` for your environment.

## Data sources

- **Mouse pancreas MALDI**: Zenodo <https://doi.org/10.5281/zenodo.3607915>.
- **Atherosclerosis regression MALDI** (macrophage metabolomics): mean spectra
  on Zenodo (<https://zenodo.org/records/20533948>); the per-spot spectra that the notebooks read are
  available from the data owners upon reasonable request.
- **Canine sarcoma MALDI**: ProteomeXchange Consortium (PRIDE) accession
  `PXD010990` (<https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD010990>).
