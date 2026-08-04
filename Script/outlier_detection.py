import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

os.makedirs("outputs", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

# ============================================================
# 1. CHARGER LES DONNEES
# ============================================================
counts = read_csv_auto("outputs/counts_filtered_for_deseq.csv")  # samples x genes
library_sizes = counts.sum(axis=1)
cpm = counts.div(library_sizes, axis=0) * 1e6
log2_cpm = np.log2(cpm + 1)

meta_rows = []
for sample in counts.index:
    parts = sample.split("_")
    meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1]})
metadata = pd.DataFrame(meta_rows).set_index("sample").loc[counts.index]

# ============================================================
# 2. CLUSTERING HIERARCHIQUE (dendrogramme)
# ============================================================
corr_matrix = log2_cpm.T.corr()
dist_matrix = 1 - corr_matrix
condensed = squareform(dist_matrix.values, checks=False)
Z = linkage(condensed, method="average")

fig, ax = plt.subplots(figsize=(10, 6))
dendrogram(Z, labels=log2_cpm.index.tolist(), ax=ax, leaf_rotation=90)
ax.set_title("Clustering hierarchique des echantillons (1 - correlation)")
ax.set_ylabel("Distance")
plt.tight_layout()
plt.savefig("outputs/outlier_dendrogram.pdf")
plt.close()
print("Dendrogramme sauvegarde: outputs/outlier_dendrogram.png")
print("-> Verifiez visuellement si un echantillon ne se regroupe PAS avec ses 4 replicats.")

# ============================================================
# 3. DISTANCE AU CENTRE DE GROUPE (sur le PCA deja calcule)
# ============================================================
try:
    pca_coords = pd.read_csv("outputs/qc_pca_coordinates.csv", index_col=0)
    print("\n--- Distance de chaque echantillon au centre (centroide) de son groupe (PC1/PC2) ---")

    results = []
    for (sex, genotype), grp in pca_coords.groupby(["sex", "genotype"]):
        centroid = grp[["PC1", "PC2"]].mean()
        dists = np.sqrt(((grp[["PC1", "PC2"]] - centroid) ** 2).sum(axis=1))
        group_mean = dists.mean()
        group_std = dists.std()
        for sample, d in dists.items():
            z = (d - group_mean) / group_std if group_std > 0 else 0
            results.append({
                "sample": sample, "sex": sex, "genotype": genotype,
                "distance_au_centroide": d, "z_score_distance": z,
                "outlier_suspect": abs(z) > 2
            })

    outlier_df = pd.DataFrame(results).sort_values("z_score_distance", ascending=False)
    outlier_df.to_csv("outputs/outlier_pca_distances.csv", index=False)
    print(outlier_df.to_string(index=False))

    suspects = outlier_df[outlier_df["outlier_suspect"]]
    if len(suspects) > 0:
        print(f"\nATTENTION: {len(suspects)} echantillon(s) suspect(s) (>2 ecarts-types du centre de leur groupe):")
        print(suspects[["sample", "z_score_distance"]].to_string(index=False))
    else:
        print("\nAucun echantillon n'est a plus de 2 ecarts-types du centre de son groupe (PCA).")
except FileNotFoundError:
    print("\nFichier outputs/qc_pca_coordinates.csv introuvable -- relancez analyse.py d'abord.")

print("\n=== TERMINE ===")
