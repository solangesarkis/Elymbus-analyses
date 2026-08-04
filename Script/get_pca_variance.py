import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

counts = pd.read_csv("outputs_full/counts_filtered_for_deseq.csv", index_col=0)
library_sizes = counts.sum(axis=1)
cpm = counts.div(library_sizes, axis=0) * 1e6
log2_cpm = np.log2(cpm + 1)
gene_var = log2_cpm.var(axis=0)
top_genes = gene_var.sort_values(ascending=False).head(500).index
pca_input = log2_cpm[top_genes]
pca = PCA(n_components=10)
pca.fit(pca_input.values)
explained = pca.explained_variance_ratio_ * 100
print(f"PC1: {explained[0]:.1f}%")
print(f"PC2: {explained[1]:.1f}%")