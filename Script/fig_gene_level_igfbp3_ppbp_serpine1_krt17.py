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

os.makedirs("outputs_full/figures", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

# ============================================================
# 1. CHARGER LES COMPTAGES ET LES GENES D'INTERET
# ============================================================
counts_for_deseq = pd.read_csv("outputs_full/counts_filtered_for_deseq.csv", index_col=0)
counts_genes_x_samples = counts_for_deseq.T
library_sizes = counts_genes_x_samples.sum(axis=0)
cpm = counts_genes_x_samples.div(library_sizes, axis=1) * 1e6
log2_cpm = np.log2(cpm + 1)

ref_file = "data/Raw count matrices Cornea/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ")

genes_of_interest = ["Igfbp3", "Ppbp", "Serpine1", "Krt17"]
gene_ids = gene_symbols[gene_symbols["Gene Symbol"].isin(genes_of_interest)]

meta_rows = []
for sample in log2_cpm.columns:
    parts = sample.split("_")
    meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1], "treatment": parts[2]})
metadata = pd.DataFrame(meta_rows).set_index("sample")

# ============================================================
# 2. ORDRE DES GROUPES SUR L'AXE X
# ============================================================
group_order = []
for sex in ["Female", "Male"]:
    for genotype in ["Healthy", "Mutant"]:
        for treatment in ["Untreated", "BIM+BAK", "BIM-PF"]:
            group_order.append((sex, genotype, treatment))

group_labels = []
for sex, genotype, treatment in group_order:
    geno_label = "Pitx2$^{+/+}$" if genotype == "Healthy" else "Pitx2$^{egl1/egl1}$"
    group_labels.append(f"{sex[0]} {geno_label} {treatment}")

treatment_colors = {"Untreated": "#999999", "BIM+BAK": "#2166ac", "BIM-PF": "#d62728"}

# ============================================================
# 3. GRAPHIQUE : dimensions ajustees pour A4, moitie inferieure
#    A4 = 21 x 29.7 cm -> moitie = ~21 x 14.85 cm -> avec marges
#    utiles ~ 19 x 12 cm = 7.5 x 4.7 pouces
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(7.5, 4.2))

for ax, gene_symbol in zip(axes, genes_of_interest):
    matching = gene_ids[gene_ids["Gene Symbol"] == gene_symbol]
    if len(matching) == 0:
        ax.set_title(f"{gene_symbol}\n(not found)", fontsize=7)
        continue
    gene_id = matching.index[0]
    expr = log2_cpm.loc[gene_id]

    for i, (sex, genotype, treatment) in enumerate(group_order):
        samples = metadata[(metadata["sex"] == sex) & (metadata["genotype"] == genotype) &
                            (metadata["treatment"] == treatment)].index
        vals = expr[samples]
        jitter = np.random.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter([i + j for j in jitter], vals, color=treatment_colors[treatment],
                   s=14, zorder=3, edgecolors="black", linewidths=0.3)
        ax.scatter([i], [vals.mean()], color=treatment_colors[treatment], marker="_",
                   s=140, zorder=4, linewidths=1.4)

    for sep in [2.5, 5.5, 8.5]:
        ax.axvline(sep, color="gray", linewidth=0.4, linestyle=":")

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_labels, fontsize=3.8, rotation=90)
    ax.set_ylabel("log$_2$(CPM+1)", fontsize=6.5)
    ax.tick_params(axis="y", labelsize=5.5)
    ax.set_title(gene_symbol, fontsize=8, fontweight="bold", style="italic")

from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markeredgecolor='black',
                           markersize=4, label=t) for t, c in treatment_colors.items()]
fig.legend(handles=legend_elements, loc="upper center", ncol=3, fontsize=6,
           bbox_to_anchor=(0.5, 1.08), frameon=False)

fig.suptitle("Animal-level normalized counts for Igfbp3, Ppbp, Serpine1, and Krt17",
              fontsize=8, y=1.16)

plt.tight_layout()
plt.savefig("outputs_full/figures/Figure_gene_level_interaction_genes_A4half.pdf",
            bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_full/figures/Figure_gene_level_interaction_genes_A4half.pdf")