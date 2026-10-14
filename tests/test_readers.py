"""Readers for SCiLS CSV exports and imzML files, on the bundled toy slice."""

import numpy as np
import pandas as pd
import pytest

import pgalign as pg


@pytest.fixture(scope="module")
def toy():
    return pg.datasets.load_mouse_pancreas_toy()


def test_read_scils_csv_round_trip(toy, tmp_path):
    unk, ref = toy
    mz = unk.var["m/z"].to_numpy(dtype=float)
    X = np.asarray(unk.X, dtype=float)  # a SCiLS export is text; float32 would not survive the round trip
    spectra = pd.DataFrame(X, columns=[str(v) for v in mz])
    spectra.insert(0, "m/z", [f"Spot {i + 1}" for i in range(unk.n_obs)])
    spectra.to_csv(tmp_path / "spectra.csv", sep=";", index=False)
    spots = pd.DataFrame({"Spot index": range(unk.n_obs),
                          "x": unk.obs["x"].to_numpy(dtype=float) + 0.5,
                          "y": unk.obs["y"].to_numpy(dtype=float)})
    spots.to_csv(tmp_path / "spots.csv", sep=";", index=False)

    adata = pg.read_scils_csv(tmp_path / "spectra.csv", tmp_path / "spots.csv")
    # pandas 3 parses decimal text to within a unit in the last place, not always exactly
    np.testing.assert_allclose(np.asarray(adata.X), X, rtol=1e-13, atol=0)
    np.testing.assert_array_equal(adata.var["m/z"].to_numpy(), mz)
    assert adata.obs_names[0] == "Spot 1"
    np.testing.assert_array_equal(adata.obs["x"].to_numpy(), spots["x"].to_numpy())
    assert list(pg.align(adata, ref).shifts["offset"]) == list(pg.align(unk, ref).shifts["offset"])


def test_read_scils_csv_rejects_mismatched_spots(toy, tmp_path):
    unk, _ = toy
    spectra = pd.DataFrame(np.asarray(unk.X)[:3], columns=[str(v) for v in unk.var["m/z"]])
    spectra.insert(0, "m/z", ["Spot 1", "Spot 2", "Spot 3"])
    spectra.to_csv(tmp_path / "spectra.csv", sep=";", index=False)
    pd.DataFrame({"Spot index": [0, 1], "x": [0.0, 1.0], "y": [0.0, 0.0]}).to_csv(
        tmp_path / "spots.csv", sep=";", index=False)
    with pytest.raises(ValueError, match="2 spots but"):
        pg.read_scils_csv(tmp_path / "spectra.csv", tmp_path / "spots.csv")


def _write_imzml(path, adata, mode, drop_last_bin_of_first=False):
    from pyimzml.ImzMLWriter import ImzMLWriter

    mz = adata.var["m/z"].to_numpy(dtype=float)
    X = np.asarray(adata.X, dtype=float)
    with ImzMLWriter(str(path), mode=mode, intensity_dtype=np.float64) as writer:
        for i in range(adata.n_obs):
            keep = slice(0, -1) if (drop_last_bin_of_first and i == 0) else slice(None)
            writer.addSpectrum(mz[keep], X[i, keep], (int(adata.obs["x"].iloc[i]), int(adata.obs["y"].iloc[i]), 1))


def test_read_imzml_continuous(toy, tmp_path):
    pytest.importorskip("pyimzml")
    unk, _ = toy
    _write_imzml(tmp_path / "toy.imzML", unk, "continuous")
    adata = pg.read_imzml(tmp_path / "toy.imzML")
    np.testing.assert_array_equal(np.asarray(adata.X), np.asarray(unk.X))
    np.testing.assert_array_equal(adata.var["m/z"].to_numpy(), unk.var["m/z"].to_numpy(dtype=float))
    np.testing.assert_array_equal(adata.obs["x"].to_numpy(), unk.obs["x"].to_numpy())


def test_read_imzml_processed_needs_and_uses_a_grid(toy, tmp_path):
    pytest.importorskip("pyimzml")
    unk, _ = toy
    _write_imzml(tmp_path / "toy.imzML", unk, "processed", drop_last_bin_of_first=True)
    with pytest.raises(ValueError, match="pass increment="):
        pg.read_imzml(tmp_path / "toy.imzML")
    adata = pg.read_imzml(tmp_path / "toy.imzML", increment=0.25, mz_range=(100, 250.25))
    np.testing.assert_array_equal(adata.var["m/z"].to_numpy(), unk.var["m/z"].to_numpy(dtype=float))
    X = np.asarray(unk.X)
    # each grid point carries a recorded value, so gridding returns it; only the dropped bin is filled in
    np.testing.assert_allclose(np.asarray(adata.X)[1:], X[1:], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(np.asarray(adata.X)[0, :-1], X[0, :-1], rtol=1e-12, atol=1e-12)
