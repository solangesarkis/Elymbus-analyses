import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

def relabel(sample_name):
    label = sample_name.replace("Healthy", "Pitx2$^{+/+}$")
    label = label.replace("Mutant", "Pitx2$^{egl1/egl1}$")
    return label

os.makedirs("outputs_TG_final/figures", exist_ok=True)

# ============================================================
# 1. CHARGER LES DONNEES FINALES (50 echantillons, post-nettoyage)
# ============================================================
counts_for_deseq = read_csv_auto("outputs_TG_final/counts_filtered_for_deseq.csv")  # samples x genes
counts_rounded = counts_for_deseq.T  # genes x samples

meta_rows = []
for sample in counts_rounded.columns:
    parts = sample.split("_")
    meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1], "treatment": parts[2]})
metadata = pd.DataFrame(meta_rows).set_index("sample")

library_sizes = counts_rounded.sum(axis=0)
cpm = counts_rounded.div(library_sizes.replace(0, np.nan), axis=1) * 1e6
cpm = cpm.fillna(0)
log2_cpm = np.log2(cpm + 1)

# ============================================================
# 2. HEATMAP DE CORRELATION
# ============================================================
corr_matrix = log2_cpm.corr()
off_diag_mask = ~np.eye(corr_matrix.shape[0], dtype=bool)
vmin_auto = corr_matrix.values[off_diag_mask].min()
print(f"Correlation min observee (hors diagonale): {vmin_auto:.4f}")

display_labels = [relabel(s) for s in corr_matrix.columns]

fig, ax = plt.subplots(figsize=(12, 11))
im = ax.imshow(corr_matrix.values, vmin=vmin_auto, vmax=1.0, cmap="viridis")
ax.set_xticks(range(len(corr_matrix.columns)))
ax.set_xticklabels(display_labels, rotation=90, fontsize=5)
ax.set_yticks(range(len(corr_matrix.columns)))
ax.set_yticklabels(display_labels, fontsize=5)
plt.colorbar(im, ax=ax, label="Correlation de Pearson")
ax.set_title("Correlation entre echantillons - TG traitement (50 ech., log2 CPM+1)")
plt.tight_layout()
plt.savefig("outputs_TG_final/figures/qc_correlation_heatmap_TG_treatment.pdf")
plt.close()
print("Sauvegarde: outputs_TG_final/figures/qc_correlation_heatmap_TG_treatment.pdf")

# ============================================================
# 3. PCA
# ============================================================
gene_var = log2_cpm.var(axis=1)
top_genes = gene_var.sort_values(ascending=False).head(500).index
pca_input = log2_cpm.loc[top_genes].T
pca = PCA(n_components=10)
pca_coords = pca.fit_transform(pca_input.values)
explained_var = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(pca_coords[:, :2], columns=["PC1", "PC2"], index=pca_input.index)
pca_df = pca_df.join(metadata)
pca_df.to_csv("outputs_TG_final/qc_pca_coordinates_treatment.csv")
print(f"\nPC1: {explained_var[0]:.1f}%, PC2: {explained_var[1]:.1f}%")

fig, ax = plt.subplots(figsize=(10, 8))
colors = {"Healthy": "tab:blue", "Mutant": "tab:red"}
markers = {"Female": "o", "Male": "^"}
sizes = {"Untreated": 40, "BIM+BAK": 90, "BIM-PF": 160}

genotype_labels = {"Healthy": "Pitx2$^{+/+}$", "Mutant": "Pitx2$^{egl1/egl1}$"}

for (sex, genotype, treatment), grp in pca_df.groupby(["sex", "genotype", "treatment"]):
    label = f"{sex}-{genotype_labels[genotype]}-{treatment}"
    ax.scatter(grp["PC1"], grp["PC2"], label=label,
               c=colors[genotype], marker=markers[sex], s=sizes[treatment], alpha=0.7,
               edgecolors="black", linewidths=0.5)

ax.set_xlabel(f"PC1 ({explained_var[0]:.1f}% variance)")
ax.set_ylabel(f"PC2 ({explained_var[1]:.1f}% variance)")
ax.set_title("PCA - TG treatment (n=50, post-decontamination)")
ax.legend(fontsize=6, bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("outputs_TG_final/figures/qc_pca_TG_treatment.pdf")
plt.close()
print("Sauvegarde: outputs_TG_final/figures/qc_pca_TG_treatment.pdf")

print("\n=== TERMINE ===")
