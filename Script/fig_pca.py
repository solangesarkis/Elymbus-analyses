import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

os.makedirs("outputs_full/figures", exist_ok=True)

pca_df = pd.read_csv("outputs_full/qc_pca_coordinates.csv", index_col=0)

# Recuperer le % de variance depuis le fichier existant si possible, sinon valeurs connues
pc1_var, pc2_var = 38.8, 13.2

def confidence_ellipse(x, y, ax, n_std=1.5, **kwargs):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                       facecolor="none", **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)
    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)

treatments = ["Untreated", "BIM+BAK", "BIM-PF"]
colors = {"Healthy": "#1f77b4", "Mutant": "#d62728"}
markers = {"Female": "o", "Male": "^"}

fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharex=True, sharey=True)

for ax, treatment in zip(axes, treatments):
    sub = pca_df[pca_df["treatment"] == treatment]
    for (sex, genotype), grp in sub.groupby(["sex", "genotype"]):
        ax.scatter(grp["PC1"], grp["PC2"],
                   c=colors[genotype], marker=markers[sex], s=90,
                   edgecolors="black", linewidths=0.6, alpha=0.85,
                   label=f"{sex} {genotype}")
        confidence_ellipse(grp["PC1"].values, grp["PC2"].values, ax,
                            edgecolor=colors[genotype], linestyle="--", linewidth=1, alpha=0.5)
    ax.set_title(treatment, fontsize=12)
    ax.set_xlabel(f"PC1 ({pc1_var:.1f}%)")
    ax.axhline(0, color="gray", linewidth=0.4)
    ax.axvline(0, color="gray", linewidth=0.4)

axes[0].set_ylabel(f"PC2 ({pc2_var:.1f}%)")

# Legende unique, partagee, en dehors des panneaux
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.08),
           ncol=4, frameon=False, fontsize=10)

fig.suptitle("PCA by treatment (color = genotype, shape = sex, ellipse = group spread)",
             fontsize=11, y=1.15)

plt.tight_layout()
plt.savefig("outputs_full/figures/fig_pca_improved.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_full/figures/fig_pca_improved.pdf")
