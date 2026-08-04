import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

os.makedirs("outputs_full/figures", exist_ok=True)

comp = pd.read_csv("outputs_full/gsea_interaction_comparison_BAK_vs_PF.csv", sep=None, engine="python")
shared = comp.dropna(subset=["NES_BIMplusBAK", "NES_BIMminusPF"]).copy()
print(f"Voies partagees: {len(shared)}")

shared["Term_short"] = shared["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True)
shared = shared.sort_values("NES_BIMplusBAK")

fig, ax = plt.subplots(figsize=(7, max(4, len(shared) * 0.5)))

for col_nes, col_fdr, xpos, label in [
    ("NES_BIMplusBAK", "FDR_BIMplusBAK", 0, "BIM+BAK"),
    ("NES_BIMminusPF", "FDR_BIMminusPF", 1, "BIM-PF"),
]:
    sizes = -np.log10(shared[col_fdr].clip(lower=1e-10)) * 60
    sc = ax.scatter([xpos] * len(shared), range(len(shared)), s=sizes,
                     c=shared[col_nes], cmap="RdBu_r", vmin=-3, vmax=3,
                     edgecolors="black", linewidths=0.6)

ax.set_yticks(range(len(shared)))
ax.set_yticklabels(shared["Term_short"], fontsize=9)
ax.set_xticks([0, 1])
ax.set_xticklabels(["BIM+BAK", "BIM-PF"])
ax.set_xlim(-0.5, 1.5)
ax.set_title(f"Pathways significant in BOTH treatments (n={len(shared)})\n(size = -log10 FDR, color = NES)", fontsize=10)

cbar = plt.colorbar(sc, ax=ax, shrink=0.5)
cbar.set_label("NES")

plt.tight_layout()
plt.savefig("outputs_full/figures/fig4b_gsea_shared_pathways.pdf")
plt.close()
print("Sauvegarde: outputs_full/figures/fig4b_gsea_shared_pathways.pdf")
