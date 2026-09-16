# Bicomposite Nanotherapeutic Architecture — v7.2 candidate

Audited computational research package derived from Zenodo v7.1 (DOI: 10.5281/zenodo.22035446).

This candidate separates three evidence layers: **public human-source transcriptomics** (GEO GSE274103; expected sample input GSM8443450), **derived bioinformatics output** (`cluster_7_differential_markers.csv`), and **synthetic/in-silico pharmacodynamic data** (`PDAC_In_Silico_Simulation_Data(101).xlsx`).

The supplied v7.1 archive does not contain experimental provenance establishing the 100 workbook rows as measurements from 100 individual patient-derived organoids. They are therefore described here as simulated/in-silico instances. The deposited top-20 Cluster 7 marker table does not contain MMP9, so this candidate does not automatically label Cluster 7 an “MMP9 invasive niche”.

## Scope
Computational proof of concept and reproducibility scaffold. It is not clinical validation, a clinical decision tool, or evidence of therapeutic efficacy or safety.

## Reproduction
Install `requirements.txt`, obtain the required GEO input separately, then run:

    python src/analyze_cluster_7.py --matrix /path/to/GSM8443450_PDAC-p2_filtered_feature_bc_matrix.h5

See `docs/` for provenance, methods, audit findings and limitations.
