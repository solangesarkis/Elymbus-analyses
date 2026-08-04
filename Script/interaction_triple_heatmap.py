import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

os.makedirs("outputs_full/figures", exist_ok=True)

df = pd.read_csv("outputs_full/interaction_triple_BIMPF_annotated.csv", sep=None, engine="python")
df = df.sort_values("padj")

lfc_cols = ["LFC_Female_Healthy", "LFC_Female_Mutant", "LFC_Male_Healthy", "LFC_Male_Mutant"]
labels = ["Female\nHealthy", "Female\nMutant", "Male\nHealthy", "Male\nMutant"]

top_n = min(40, len(df))
mat = df[lfc_cols].head(top_n).values
gene_labels = df["Gene Symbol"].head(top_n).tolist()

fig, ax = plt.subplots(figsize=(6, 12))
vmax = np.nanmax(np.abs(mat))
im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(4))
ax.set_xticklabels(labels)
ax.set_yticks(range(top_n))
ax.set_yticklabels(gene_labels, fontsize=7)
ax.set_title(f"Top {top_n} triple-interaction genes\n(BIM-PF vs Untreated, LFC by subgroup)", fontsize=10)
cbar = plt.colorbar(im, ax=ax, shrink=0.5)
cbar.set_label("log2FoldChange\n(BIM-PF vs Untreated)")
plt.tight_layout()
plt.savefig("outputs_full/figures/interaction_triple_heatmap.pdf")
plt.close()
print(f"Sauvegarde: outputs_full/figures/interaction_triple_heatmap.pdf ({top_n} genes affiches)")
