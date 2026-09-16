# Reproducibility

The supplied v7.1 script used a Windows-specific absolute path and automatic runtime package installation. This candidate removes both. The corrected script receives the matrix with `--matrix` and output location with `--outdir`. The required GEO HDF5 input was not included in the six-file Zenodo v7.1 archive supplied for audit, so full end-to-end numerical regeneration of the deposited marker table and UMAP cannot yet be claimed from that archive alone.
