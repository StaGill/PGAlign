"""Command line interface: ``pgalign align UNKNOWN REFERENCE -o ALIGNED.h5ad``."""

import argparse
import sys
from pathlib import Path

from . import __version__
from .pipeline import align
from .readers import read_h5ad, read_imzml, read_scils_csv, write_h5ad


def _read(path, spots, increment, mz_range):
    suffix = Path(path).suffix.lower()
    if suffix == ".imzml":
        return read_imzml(path, increment=increment, mz_range=mz_range)
    if suffix == ".h5ad":
        return read_h5ad(path)
    if suffix == ".csv":
        return read_scils_csv(path, spots)
    raise ValueError(f"cannot read {path}: expected an .imzML, .h5ad or SCiLS .csv file")


def _parser():
    parser = argparse.ArgumentParser(
        prog="pgalign",
        description="PGAlign: peak group alignment of mass spectrometry data.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    cmd = commands.add_parser(
        "align",
        help="align an unknown dataset onto a reference",
        description=(
            "Align UNKNOWN onto REFERENCE with the full PGAlign pipeline. Each input is an "
            ".imzML, .h5ad or SCiLS Lab .csv export. Writes the aligned unknown spectra to "
            "OUT, and next to it the reference on the same bins (<stem>_reference.h5ad), "
            "the per-group shifts (<stem>_shifts.csv) and the ambiguity table "
            "(<stem>_ambiguity.csv)."
        ),
    )
    cmd.add_argument("unknown", help="the dataset to align")
    cmd.add_argument("reference", help="the dataset to align onto")
    cmd.add_argument("-o", "--out", required=True, help="path of the aligned .h5ad to write")
    cmd.add_argument("--unknown-spots", help="Region Spots.csv for a SCiLS .csv unknown")
    cmd.add_argument("--reference-spots", help="Region Spots.csv for a SCiLS .csv reference")
    cmd.add_argument("--increment", type=float,
                     help="grid spacing in m/z for imzML input; required for processed imzML")
    cmd.add_argument("--mz-range", type=float, nargs=2, metavar=("START", "END"),
                     help="m/z range of that grid, END excluded")
    cmd.add_argument("--no-smooth", dest="smooth", action="store_false",
                     help="call peaks with PeakCalling instead of the moving average")
    cmd.add_argument("--alpha", type=float, default=0.9, help="peak calling quantile (default 0.9)")
    cmd.add_argument("--percentile", type=float, default=0.9,
                     help="peak grouping parameter (default 0.9)")
    cmd.add_argument("--delta", type=float, default=0, help="pairing criterion (default 0)")
    cmd.add_argument("--tau", type=float, default=0.1,
                     help="tolerance of the ambiguity diagnostic (default 0.1)")
    cmd.add_argument("--no-normalize", dest="normalize", action="store_false",
                     help="skip the per-spectrum normalization")
    return parser


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        unk = _read(args.unknown, args.unknown_spots, args.increment, args.mz_range)
        ref = _read(args.reference, args.reference_spots, args.increment, args.mz_range)
        res = align(unk, ref, smooth=args.smooth, alpha=args.alpha, percentile=args.percentile,
                    delta=args.delta, tau=args.tau, normalize=args.normalize)
    except (ValueError, ImportError) as err:
        parser.exit(2, f"pgalign: error: {err}\n")
    out = Path(args.out)
    stem = out.with_suffix("")
    write_h5ad(res.aligned, out)
    write_h5ad(res.reference, f"{stem}_reference.h5ad")
    res.shifts.to_csv(f"{stem}_shifts.csv", index=False)
    res.ambiguity.to_csv(f"{stem}_ambiguity.csv", index=False)
    print(res.summary())
    return 0


if __name__ == "__main__":
    sys.exit(main())
