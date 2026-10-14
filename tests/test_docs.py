"""Run the Python examples of README.md and docs/ as written, on the bundled toy data."""

import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PAGES = ["README.md"] + sorted(str(p.relative_to(ROOT)) for p in (ROOT / "docs").glob("*.md"))
# Examples that read the reader's own files cannot run here.
NEEDS_FILES = ('read_imzml("', 'read_scils_csv("', 'read_h5ad("')
# The README step-by-step example leaves the inputs to the reader; supply the toy data.
PLACEHOLDER = ("MALDIdataAnn1 = ...   # unknown spectrum\n", "MALDIdataAnn2 = ...   # reference spectrum\n")
TOY = "MALDIdataAnn1, MALDIdataAnn2 = pg.datasets.load_mouse_pancreas_toy()\n"


def _blocks(page):
    text = (ROOT / page).read_text(encoding="utf-8")
    return [b for b in re.findall(r"```python\n(.*?)```", text, flags=re.S) if not any(s in b for s in NEEDS_FILES)]


def _run(page):
    namespace = {}
    for block in _blocks(page):
        exec(block.replace(PLACEHOLDER[0], TOY).replace(PLACEHOLDER[1], ""), namespace)
        if "res" in namespace and "res.summary()" in block:
            assert namespace["res"].summary() in block, f"{page}: the summary in a comment no longer matches"
    return namespace


@pytest.mark.parametrize("page", [p for p in PAGES if (ROOT / p).exists() and _blocks(p)])
def test_examples_run_as_written(page):
    _run(page)


def test_readme_examples_align_onto_the_reference_grid():
    if not (ROOT / "README.md").exists():
        pytest.skip("README.md is not in this checkout")
    ns = _run("README.md")
    aligned, ref = ns["MALDIdataAnn1_aligned"], ns["MALDIdataAnn2"]
    np.testing.assert_array_equal(aligned.var["m/z"].to_numpy(), ref.var["m/z"].to_numpy()[:aligned.n_vars])
    assert len(ns["table"]) > 0


def test_tutorial_step_by_step_equals_align():
    if not (ROOT / "docs" / "tutorial.md").exists():
        pytest.skip("docs/ is not in this checkout")
    ns = _run("docs/tutorial.md")
    # the toy shift is a rigid +2, so the matrix alone cannot tell the two peak callers apart
    np.testing.assert_array_equal(ns["align"].align_group, ns["res"].steps.align_group)
    np.testing.assert_array_equal(np.asarray(ns["aligned"].X), np.asarray(ns["res"].aligned.X))
