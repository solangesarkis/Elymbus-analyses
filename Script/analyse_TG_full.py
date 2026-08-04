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
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

os.makedirs("outputs_TG_full", exist_ok=True)

# ============================================================
# 1. CHARGEMENT DES 12 FICHIERS
# ============================================================
files = {
    ("Female", "Healthy", "Untreated"): "data/Female Control NT TG.xlsx",
    ("Female", "Mutant",  "Untreated"): "data/Female Mutant NT TG.xlsx",
    ("Male",   "Healthy", "Untreated"): "data/Male Control NT TG.xlsx",
    ("Male",   "Mutant",  "Untreated"): "data/Male Mutant NT TG.xlsx",
    ("Female", "Healthy", "BIM+BAK"):   "data/Female Control BIM+BAK TG.xlsx",
    ("Female", "Mutant",  "BIM+BAK"):   "data/Female Mutant BIM+BAK TG.xlsx",
    ("Male",   "Healthy", "BIM+BAK"):   "data/Male Control BIM+BAK TG.xlsx",
    ("Male",   "Mutant",  "BIM+BAK"):   "data/Male Mutant BIM+BAK TG.xlsx",
    ("Female", "Healthy", "BIM-PF"):    "data/Female Control BIM-PF TG.xlsx",
    ("Female", "Mutant",  "BIM-PF"):    "data/Female Mutant BIM-PF TG.xlsx",
    ("Male",   "Healthy", "BIM-PF"):    "data/Male Control BIM-PF TG.xlsx",
    ("Male",   "Mutant",  "BIM-PF"):    "data/Male Mutant BIM-PF TG.xlsx",
}

all_counts_raw = []
all_gene_info = []
metadata_rows = []

for (sex, genotype, treatment), path in files.items():
    df = pd.read_excel(path)
    df = df.set_index("Gene ID")
    id_cols = [c for c in ("Gene Symbol", "Type") if c in df.columns]
    all_gene_info.append(df[id_cols])
    count_cols = [c for c in df.columns if c not in id_cols]
    counts_only = df[count_cols]
    new_names = [f"{sex}_{genotype}_{treatment}_{i+1}" for i in range(counts_only.shape[1])]
    counts_only.columns = new_names
    all_counts_raw.append(counts_only)
    for name in new_names:
        metadata_rows.append({"sample": name, "sex": sex, "genotype": genotype, "treatment": treatment})

counts_raw = pd.concat(all_counts_raw, axis=1, join="outer").fillna(0)
gene_info = pd.concat(all_gene_info).groupby(level=0).first().reindex(counts_raw.index)

metadata = pd.DataFrame(metadata_rows).set_index("sample").loc[counts_raw.columns]
print("--- Metadonnees (12 groupes x 5 replicats = 60 echantillons) ---")
print(metadata)
print(f"\nTotal genes uniques: {counts_raw.shape[0]}")

# ============================================================
# 2. ARRONDI (0.5 vers le haut) + SAUVEGARDE AUDIT
# ============================================================
counts_raw.to_csv("outputs_TG_full/counts_raw_unrounded.csv")
counts_rounded = np.floor(counts_raw.values + 0.5)
counts_rounded = pd.DataFrame(counts_rounded, index=counts_raw.index, columns=counts_raw.columns).astype(int)
counts_rounded.to_csv("outputs_TG_full/counts_rounded.csv")

# ============================================================
# 3. FILTRAGE : >=10 reads dans >=5 echantillons (sur les 60)
# ============================================================
keep = (counts_rounded >= 10).sum(axis=1) >= 5
counts_filtered = counts_rounded.loc[keep]
print(f"\nGenes avant filtrage: {counts_rounded.shape[0]}")
print(f"Genes apres filtrage: {counts_filtered.shape[0]}")
counts_for_deseq = counts_filtered.T
counts_for_deseq.to_csv("outputs_TG_full/counts_filtered_for_deseq.csv")

# ============================================================
# 4. QC : tailles de bibliotheque, correlation, PCA
# ============================================================
library_sizes = counts_rounded.sum(axis=0)
print("\n--- Tailles de bibliotheque (min/max/mean) ---")
print(f"Min: {library_sizes.min():.0f}, Max: {library_sizes.max():.0f}, Mean: {library_sizes.mean():.0f}")

fig, ax = plt.subplots(figsize=(14, 6))
library_sizes.plot(kind="bar", ax=ax)
ax.set_ylabel("Total de lectures")
ax.set_title("Tailles de bibliotheque - TG - 60 echantillons")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()
plt.savefig("outputs_TG_full/qc_library_sizes.pdf")
plt.close()

cpm = counts_rounded.div(library_sizes.replace(0, np.nan), axis=1) * 1e6
cpm = cpm.fillna(0)
log2_cpm = np.log2(cpm + 1)

corr_matrix = log2_cpm.corr()
fig, ax = plt.subplots(figsize=(12, 11))
im = ax.imshow(corr_matrix.values, vmin=0, vmax=1.0, cmap="viridis")
ax.set_xticks(range(len(corr_matrix.columns)))
ax.set_xticklabels(corr_matrix.columns, rotation=90, fontsize=5)
ax.set_yticks(range(len(corr_matrix.columns)))
ax.set_yticklabels(corr_matrix.columns, fontsize=5)
plt.colorbar(im, ax=ax, label="Correlation de Pearson")
ax.set_title("Correlation entre les 60 echantillons TG (log2 CPM+1)")
plt.tight_layout()
plt.savefig("outputs_TG_full/qc_sample_correlation_heatmap.pdf")
plt.close()

dist_matrix = 1 - corr_matrix
condensed = squareform(dist_matrix.values, checks=False)
Z = linkage(condensed, method="average")
fig, ax = plt.subplots(figsize=(16, 6))
dendrogram(Z, labels=log2_cpm.columns.tolist(), ax=ax, leaf_rotation=90, leaf_font_size=6)
ax.set_title("Clustering hierarchique - TG - 60 echantillons")
plt.tight_layout()
plt.savefig("outputs_TG_full/qc_dendrogram.pdf")
plt.close()

gene_var = log2_cpm.var(axis=1)
top_genes = gene_var.sort_values(ascending=False).head(500).index
pca_input = log2_cpm.loc[top_genes].T
pca = PCA(n_components=10)
pca_coords = pca.fit_transform(pca_input.values)
explained_var = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(pca_coords[:, :2], columns=["PC1", "PC2"], index=pca_input.index)
pca_df = pca_df.join(metadata)
pca_df.to_csv("outputs_TG_full/qc_pca_coordinates.csv")

fig, ax = plt.subplots(figsize=(10, 8))
colors = {"Healthy": "tab:blue", "Mutant": "tab:red"}
markers = {"Female": "o", "Male": "^"}
sizes = {"Untreated": 40, "BIM+BAK": 90, "BIM-PF": 160}
for (sex, genotype, treatment), grp in pca_df.groupby(["sex", "genotype", "treatment"]):
    ax.scatter(grp["PC1"], grp["PC2"],
               label=f"{sex}-{genotype}-{treatment}",
               c=colors[genotype], marker=markers[sex], s=sizes[treatment], alpha=0.7,
               edgecolors="black", linewidths=0.5)
ax.set_xlabel(f"PC1 ({explained_var[0]:.1f}% variance)")
ax.set_ylabel(f"PC2 ({explained_var[1]:.1f}% variance)")
ax.set_title("PCA - TG - 60 echantillons")
ax.legend(fontsize=6, bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("outputs_TG_full/qc_pca.pdf")
plt.close()

# ============================================================
# 5. DETECTION D'OUTLIERS
# ============================================================
results = []
for (sex, genotype, treatment), grp in pca_df.groupby(["sex", "genotype", "treatment"]):
    centroid = grp[["PC1", "PC2"]].mean()
    dists = np.sqrt(((grp[["PC1", "PC2"]] - centroid) ** 2).sum(axis=1))
    group_mean = dists.mean()
    group_std = dists.std()
    for sample, d in dists.items():
        z = (d - group_mean) / group_std if group_std > 0 else 0
        results.append({"sample": sample, "sex": sex, "genotype": genotype, "treatment": treatment,
                         "distance_au_centroide": d, "z_score": z, "outlier_suspect": abs(z) > 2})
outlier_df = pd.DataFrame(results).sort_values("z_score", ascending=False)
outlier_df.to_csv("outputs_TG_full/outlier_detection.csv", index=False)
print("\n--- Echantillons suspects (z-score > 2) ---")
print(outlier_df[outlier_df["outlier_suspect"]].to_string(index=False))
if outlier_df["outlier_suspect"].sum() == 0:
    print("Aucun outlier detecte.")

print("\n=== QC TERMINE. Fichiers dans outputs_TG_full/ ===")
