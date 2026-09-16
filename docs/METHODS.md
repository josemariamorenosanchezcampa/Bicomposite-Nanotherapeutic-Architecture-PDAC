# Methods

The audited pipeline loads a user-supplied 10X HDF5 expression matrix, filters observations and genes, normalizes counts, applies log1p, selects highly variable genes, performs PCA, constructs a neighbor graph, computes UMAP and applies Leiden clustering. Differential expression for Cluster 7 versus the remainder uses Scanpy Wilcoxon ranking. UMAP is treated as an expression-space embedding, not an anatomical map. MMP9 is checked descriptively if present; no biological niche label is inferred automatically.
