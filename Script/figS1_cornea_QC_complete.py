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
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

os.makedirs("outputs_full/figures", exist_ok=True)

def read_csv_auto(path, **kwargs):
    try:
        return pd.read_csv(path, sep=None, engine="python", encoding="utf-8", **kwargs)
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=None, engine="python", encoding="latin-1", **kwargs)

# ============================================================
# CHARGEMENT DES DONNEES
# ============================================================
counts = pd.read_csv("outputs_full/counts_filtered_for_deseq.csv", index_col=0)
library_sizes = counts.sum(axis=1)
detected_genes = (counts > 0).sum(axis=1)
cpm = counts.div(library_sizes, axis=0) * 1e6
log2cpm = np.log2(cpm + 1)

meta_rows = []
for sample in counts.index:
    parts = sample.split("_")
    meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1], "treatment": parts[2]})
metadata = pd.DataFrame(meta_rows).set_index("sample")

GENO_DISPLAY = {"Healthy": "Pitx2+/+", "Mutant": "Pitx2egl1/egl1"}
markers = {"Female": "o", "Male": "^"}
geno_colors = {"Healthy": "#2166ac", "Mutant": "#d62728"}
treat_colors = {"Untreated": "#999999", "BIM+BAK": "#2166ac", "BIM-PF": "#d62728"}

fig, axes = plt.subplots(3, 2, figsize=(8.0, 11.0))
axA, axB, axC, axD, axE, axF = axes.flatten()

# ============================================================
# PANEL A : taille de bibliotheque
# ============================================================
lib_sorted = library_sizes.sort_values()
colors_a = [treat_colors[metadata.loc[s, "treatment"]] for s in lib_sorted.index]
axA.bar(range(len(lib_sorted)), lib_sorted.values, color=colors_a, edgecolor="black", linewidth=0.2)
axA.set_xlabel("Sample (sorted)", fontsize=6)
axA.set_ylabel("Library size", fontsize=6)
axA.tick_params(axis="both", labelsize=5)
axA.set_title("A", fontsize=9, fontweight="bold", loc="left")

# ============================================================
# PANEL B : genes detectes
# ============================================================
det_sorted = detected_genes.sort_values()
colors_b = [treat_colors[metadata.loc[s, "treatment"]] for s in det_sorted.index]
axB.bar(range(len(det_sorted)), det_sorted.values, color=colors_b, edgecolor="black", linewidth=0.2)
axB.set_xlabel("Sample (sorted)", fontsize=6)
axB.set_ylabel("Detected genes", fontsize=6)
axB.tick_params(axis="both", labelsize=5)
axB.set_title("B", fontsize=9, fontweight="bold", loc="left")

treat_legend = [Patch(facecolor=c, edgecolor="black", label=t) for t, c in treat_colors.items()]
axB.legend(handles=treat_legend, fontsize=4.5, loc="upper left")

# ============================================================
# PANEL C : PCA facettee par traitement
# ============================================================
old_pca = pd.read_csv("outputs_full/qc_pca_coordinates.csv", index_col=0)
pcs_df = old_pca[["PC1", "PC2"]].join(metadata)

PC1_VAR, PC2_VAR = 38.8, 13.2

pc1_min, pc1_max = pcs_df["PC1"].min() * 1.1, pcs_df["PC1"].max() * 1.1
pc2_min, pc2_max = pcs_df["PC2"].min() * 1.1, pcs_df["PC2"].max() * 1.1

axC.axis("off")
axC.text(-0.05, 1.08, "C", transform=axC.transAxes, fontsize=9, fontweight="bold")

for i, treat in enumerate(["Untreated", "BIM+BAK", "BIM-PF"]):
    inset = axC.inset_axes([i * 0.35, 0.05, 0.30, 0.85])
    sub = pcs_df[pcs_df["treatment"] == treat]
    for (sex, genotype), grp in sub.groupby(["sex", "genotype"]):
        inset.scatter(grp["PC1"], grp["PC2"], c=geno_colors[genotype], marker=markers[sex],
                      s=14, edgecolors="black", linewidths=0.3)
    inset.set_xlim(pc1_min, pc1_max)
    inset.set_ylim(pc2_min, pc2_max)
    inset.set_title(treat, fontsize=5.5)
    inset.set_xlabel(f"PC1 ({PC1_VAR:.1f}%)", fontsize=4.5)
    if i == 0:
        inset.set_ylabel(f"PC2 ({PC2_VAR:.1f}%)", fontsize=4.5)
    else:
        inset.set_yticklabels([])
    inset.tick_params(axis="both", labelsize=4)

geno_sex_legend = []
for sex in ["Female", "Male"]:
    for genotype in ["Healthy", "Mutant"]:
        geno_sex_legend.append(plt.Line2D([0], [0], marker=markers[sex], color="w",
                                markerfacecolor=geno_colors[genotype], markeredgecolor="black",
                                markersize=5, label=f"{sex} {GENO_DISPLAY[genotype]}"))
axC.legend(handles=geno_sex_legend, fontsize=4, loc="lower center",
           bbox_to_anchor=(0.5, -0.15), ncol=2)

# ============================================================
# PANEL D : matrice de correlation
# ============================================================
corr_matrix = log2cpm.T.corr()
im = axD.imshow(corr_matrix.values, cmap="viridis", vmin=corr_matrix.values.min(), vmax=1.0)
axD.set_xticks([]); axD.set_yticks([])
axD.set_title("D", fontsize=9, fontweight="bold", loc="left")
cbar_d = plt.colorbar(im, ax=axD, shrink=0.6)
cbar_d.set_label("Pearson r", fontsize=5.5)
cbar_d.ax.tick_params(labelsize=4.5)

# ============================================================
# PANEL E : clustering hierarchique
# ============================================================
dist_matrix = (1 - corr_matrix).copy()
dist_values = dist_matrix.values.copy()
np.fill_diagonal(dist_values, 0)
condensed = squareform(dist_values, checks=False)
Z = linkage(condensed, method="average")
dendrogram(Z, ax=axE, labels=corr_matrix.columns.tolist(), no_labels=True,
           color_threshold=0, above_threshold_color="black")
axE.set_ylabel("1 - Pearson r", fontsize=6)
axE.tick_params(axis="both", labelsize=5)
axE.set_title("E", fontsize=9, fontweight="bold", loc="left")

# ============================================================
# PANEL F : PC1-PC2 + distance au centroide
# ============================================================
pcs_df["group"] = pcs_df["sex"] + "_" + pcs_df["genotype"]
centroids = pcs_df.groupby("group")[["PC1", "PC2"]].mean()
pcs_df["dist_to_centroid"] = pcs_df.apply(
    lambda r: np.sqrt((r["PC1"] - centroids.loc[r["group"], "PC1"])**2 +
                       (r["PC2"] - centroids.loc[r["group"], "PC2"])**2), axis=1)

for (sex, genotype), grp in pcs_df.groupby(["sex", "genotype"]):
    axF.scatter(grp["PC1"], grp["PC2"], c=geno_colors[genotype], marker=markers[sex],
                s=16, edgecolors="black", linewidths=0.3, label=f"{sex} {GENO_DISPLAY[genotype]}")
    cx, cy = centroids.loc[f"{sex}_{genotype}"]
    axF.scatter(cx, cy, marker="x", s=40, color="black", linewidths=1.2, zorder=5)
axF.set_xlabel(f"PC1 ({PC1_VAR:.1f}%)", fontsize=6)
axF.set_ylabel(f"PC2 ({PC2_VAR:.1f}%)", fontsize=6)
axF.tick_params(axis="both", labelsize=5)
axF.legend(fontsize=4, loc="best")
axF.set_title("F", fontsize=9, fontweight="bold", loc="left")

plt.tight_layout()
plt.savefig("outputs_full/figures/FigureS1_cornea_QC_complete_A4.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_full/figures/FigureS1_cornea_QC_complete_A4.pdf")

max_dist = pcs_df["dist_to_centroid"].max()
print(f"\nDistance maximale au centroide du sous-groupe: {max_dist:.2f}")
print(pcs_df.sort_values("dist_to_centroid", ascending=False)[["sex", "genotype", "treatment", "dist_to_centroid"]].head(5).to_string())