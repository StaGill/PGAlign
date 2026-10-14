"""The pgalign command line on the bundled toy slice."""

import numpy as np
import pandas as pd
import pytest

import pgalign as pg
from pgalign.cli import main


@pytest.fixture
def toy_files(tmp_path):
    unk, ref = pg.datasets.load_mouse_pancreas_toy()
    pg.write_h5ad(unk, tmp_path / "unk.h5ad")
    pg.write_h5ad(ref, tmp_path / "ref.h5ad")
    return tmp_path


def test_align_writes_the_outputs(toy_files, capsys):
    out = toy_files / "aligned.h5ad"
    assert main(["align", str(toy_files / "unk.h5ad"), str(toy_files / "ref.h5ad"), "-o", str(out)]) == 0
    assert capsys.readouterr().out.strip() == "7 groups paired, 7 shifted, 0 flagged at tau = 0.1"
    aligned = pg.read_h5ad(out)
    reference = pg.read_h5ad(toy_files / "aligned_reference.h5ad")
    assert aligned.shape == reference.shape == (120, 599)
    assert set(pd.read_csv(toy_files / "aligned_shifts.csv")["offset"]) == {2}
    assert len(pd.read_csv(toy_files / "aligned_ambiguity.csv")) == 7
    res = pg.align(*pg.datasets.load_mouse_pancreas_toy())
    np.testing.assert_array_equal(np.asarray(aligned.X), np.asarray(res.aligned.X))


def test_no_smooth_uses_peakcalling(toy_files, capsys):
    main(["align", str(toy_files / "unk.h5ad"), str(toy_files / "ref.h5ad"),
          "-o", str(toy_files / "aligned.h5ad"), "--no-smooth"])
    assert capsys.readouterr().out.strip() == "15 groups paired, 15 shifted, 0 flagged at tau = 0.1"


def test_unusable_input_exits_with_the_message(toy_files, capsys):
    far = pg.read_h5ad(toy_files / "unk.h5ad")
    far.var["m/z"] = far.var["m/z"] + 1000
    pg.write_h5ad(far, toy_files / "far.h5ad")
    with pytest.raises(SystemExit) as exit_info:
        main(["align", str(toy_files / "far.h5ad"), str(toy_files / "ref.h5ad"), "-o", str(toy_files / "x.h5ad")])
    assert exit_info.value.code == 2
    assert "m/z ranges do not overlap" in capsys.readouterr().err
