# Readers

`pg.read_imzml` and `pg.read_scils_csv` turn MSI exports into the `AnnData`
objects that `pg.align` takes, and `pg.read_h5ad` reads `.h5ad` files.

::: pgalign.readers.read_imzml

::: pgalign.readers.read_scils_csv

::: pgalign.readers.read_h5ad

## Saving

Since version 0.13, anndata cannot write a column named `m/z`, so `pg.write_h5ad`
stores it as `mz` and `pg.read_h5ad` restores it.

::: pgalign.readers.write_h5ad
