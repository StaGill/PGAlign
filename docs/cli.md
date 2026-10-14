# Command line

Installing PGAlign also installs the `pgalign` command, which runs
[`pg.align`](api/pipeline.md) on two files:

```bash
pgalign align unknown.imzML reference.imzML -o aligned.h5ad --increment 0.25
```

Each input can be an `.imzML` file (read with [`pg.read_imzml`](api/readers.md),
which needs `pip install "pgalign[imzml]"`), an `.h5ad` file (read with
`pg.read_h5ad`), or a SCiLS Lab
`.csv` export (read with `pg.read_scils_csv`; pass its `Region Spots.csv` with
`--unknown-spots` or `--reference-spots`). `--increment` and `--mz-range` grid
imzML input and are required for processed imzML.

The command prints the summary line, for example
`198 groups paired, 78 shifted, 20 flagged at tau = 0.1`, and writes four files
named after `-o`; the `.h5ad` files are written with `pg.write_h5ad`, so read them
back with `pg.read_h5ad`:

| File | Content |
|---|---|
| `aligned.h5ad` | the unknown spectra on the reference m/z grid (`res.aligned`) |
| `aligned_reference.h5ad` | the reference on the same bins and scale (`res.reference`) |
| `aligned_shifts.csv` | one row per paired peak group with its offset (`res.shifts`) |
| `aligned_ambiguity.csv` | the diagnostic for ambiguous alignments (`res.ambiguity`) |

The options mirror the arguments of `pg.align`: `--no-smooth` for the
`PeakCalling` caller, `--alpha`, `--percentile`, `--delta`, `--tau`, and
`--no-normalize`. `pgalign align --help` lists them. Input that cannot be aligned
stops the command with exit status 2 and a message saying what to fix.
