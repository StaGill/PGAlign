# Reproducing the GALAXY manuscript figures

This folder contains the notebooks that generate the analyses and figures in
Deng, Zhang & Zhang (2026). Each subfolder targets one case study.

| Subfolder / file                  | Case study                                    | Manuscript section             | Dataset access                                                                 |
|-----------------------------------|-----------------------------------------------|--------------------------------|--------------------------------------------------------------------------------|
| `MousePancreate.ipynb`            | Simulation study on mouse pancreas MALDI      | Section 3.1 (Simulation study) | Zenodo: <https://doi.org/10.5281/zenodo.3607915>                               |
| `Maaike/`                         | Atherosclerosis regression (macrophages)      | Section 3.2 (Macrophage metabolomics) | Available from the data owners upon reasonable request                         |
| `Sarcomas/`                       | Canine sarcoma classification                  | Section 3.3 (Canine sarcomas)  | ProteomeXchange PRIDE accession `PXD010990`                                    |

## How to run

1. Install GALAXY from the repository root: `pip install -e .`.
2. Place the input data for the case study you want to run under a local folder
   (e.g. `data/MousePancreate/`, `data/Maaike/`, `data/CanineSarcomas/`).
3. Open the notebook; the first executable cell is a **CONFIGURE ME** block
   that sets `DATA_DIR`, `OUTPUT_DIR`, and `REPO_ROOT`. Edit those paths to
   match your environment.
4. Run the notebook end-to-end.

The notebooks were originally authored on Compute Canada (`/lustre03/...`) and
on a local macOS workstation (`/Users/anjideng1/...`). The hard-coded prefixes
have been replaced with `f"{DATA_DIR}/..."` and `f"{OUTPUT_DIR}/..."` so that
only the CONFIGURE-ME cell needs editing in a fresh environment.

The notebooks use the current manuscript's method names (`peak_calling`,
`peak_grouping`, `peak_group_pairing`, `fine_alignment_assessment`). The
legacy names (`callpeak`, `grouppeaks`, `greedy_match`, `fine_align`) remain
available as deprecated aliases and will emit `DeprecationWarning`.
