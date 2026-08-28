# Installation

GALAXY targets Python &geq; 3.9.

## From source

```bash
git clone https://github.com/StaGill/GALAXY.git
cd GALAXY
pip install -e .
```

The installed package's importable name is `GalaxyPython`:

```python
import GalaxyPython as gx
print(gx.__version__)
```

## Dependencies

GALAXY depends on the standard scientific-Python stack plus `scanpy` /
`anndata` for MALDI data containers:

- `numpy >= 1.20`
- `scipy >= 1.7`
- `pandas >= 1.3`
- `scanpy >= 1.9`
- `anndata >= 0.8`
- `tqdm`

For development (running the test suite):

```bash
pip install -e ".[test]"
pytest -q
```

## Building the docs locally

The documentation site is built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
and [mkdocstrings](https://mkdocstrings.github.io/):

```bash
pip install mkdocs-material mkdocstrings[python]
mkdocs serve
```

Then open <http://127.0.0.1:8000>.
