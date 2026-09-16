from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def run(matrix: Path, outdir: Path, resolution: float=0.8, cluster: str="7"):
    outdir.mkdir(parents=True, exist_ok=True)
    adata=sc.read_10x_h5(matrix); adata.var_names_make_unique()
    sc.pp.filter_cells(adata,min_genes=100); sc.pp.filter_genes(adata,min_cells=3)
    sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=min(2000,adata.n_vars))
    sc.pp.pca(adata,use_highly_variable=True)
    sc.pp.neighbors(adata,n_neighbors=15,n_pcs=min(30,adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata); sc.tl.leiden(adata,resolution=resolution,key_added="leiden_clusters")
    sc.pl.umap(adata,color="leiden_clusters",title="Leiden clusters (transcriptomic UMAP)",show=False)
    plt.savefig(outdir/"umap_clusters.png",dpi=300,bbox_inches="tight"); plt.close()
    available=set(map(str,adata.obs["leiden_clusters"].unique()))
    if cluster not in available: raise ValueError(f"Cluster {cluster} not found. Available: {sorted(available)}")
    sc.tl.rank_genes_groups(adata,groupby="leiden_clusters",groups=[cluster],reference="rest",method="wilcoxon")
    df=sc.get.rank_genes_groups_df(adata,group=cluster).head(20).rename(columns={"names":"gene","scores":"score","logfoldchanges":"logfoldchange"})
    cols=[c for c in ["gene","score","logfoldchange","pvals_adj"] if c in df.columns]
    df[cols].to_csv(outdir/f"cluster_{cluster}_differential_markers.csv",index=False)
    top10=df.head(10).sort_values("logfoldchange")
    plt.figure(figsize=(10,6)); plt.barh(top10["gene"],top10["logfoldchange"]); plt.xlabel(f"Log fold change (Cluster {cluster} vs rest)"); plt.ylabel("Genes"); plt.title(f"Top upregulated genes in Cluster {cluster}"); plt.tight_layout(); plt.savefig(outdir/f"cluster_{cluster}_top_genes.png",dpi=300); plt.close()
    if "MMP9" in adata.var_names:
        mask=adata.obs["leiden_clusters"].astype(str).to_numpy()==cluster
        x=adata[:,["MMP9"]].X; x=x.toarray().ravel() if hasattr(x,"toarray") else np.asarray(x).ravel()
        pd.DataFrame([{"gene":"MMP9","cluster":cluster,"mean_log_normalized_cluster":float(x[mask].mean()),"mean_log_normalized_rest":float(x[~mask].mean())}]).to_csv(outdir/"mmp9_descriptive_check.csv",index=False)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--matrix",type=Path,required=True); p.add_argument("--outdir",type=Path,default=Path("outputs/recomputed")); p.add_argument("--resolution",type=float,default=0.8); p.add_argument("--cluster",default="7"); a=p.parse_args(); run(a.matrix,a.outdir,a.resolution,a.cluster)
if __name__=="__main__": main()
