# Case studies

The notebooks under [`CodeInPaper/`](https://github.com/StaGill/GALAXY/tree/main/CodeInPaper)
reproduce the analyses in the GALAXY manuscript.

| Subfolder / file                  | Case study                              | Manuscript section             | Dataset access                                                                 |
|-----------------------------------|-----------------------------------------|--------------------------------|--------------------------------------------------------------------------------|
| `MousePancreate.ipynb`            | Simulation on mouse pancreas MALDI       | Section 3.1                    | Zenodo: <https://doi.org/10.5281/zenodo.3607915>                               |
| `Maaike/`                         | Atherosclerosis regression macrophages   | Section 3.2                    | Available from the data owners upon reasonable request                         |
| `Sarcomas/`                       | Canine sarcoma classification             | Section 3.3                    | ProteomeXchange PRIDE accession `PXD010990`                                    |

Each subfolder has its own `README.md` with dataset access details and a
**CONFIGURE ME** cell at the top of each notebook that sets `DATA_DIR`,
`OUTPUT_DIR`, and `REPO_ROOT` for your environment.

## Data sources

- **Mouse pancreas MALDI**: Zenodo <https://doi.org/10.5281/zenodo.3607915>.
- **Atherosclerosis regression MALDI** (macrophage metabolomics): available
  from the data owners upon reasonable request.
- **Canine sarcoma MALDI**: ProteomeXchange Consortium (PRIDE) accession
  `PXD010990` (<https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD010990>).
