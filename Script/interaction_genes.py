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

# ============================================================
# 1. CHARGER LES RESULTATS D'INTERACTION ET GARDER LES SIGNIFICATIFS
# ============================================================
interaction_df = pd.read_csv("outputs/DE_interaction_sex_genotype.csv", index_col=0, sep=";")
sig_genes = interaction_df[interaction_df["padj"] < 0.05].sort_values("padj")
print(f"Nombre de genes d'interaction significatifs: {len(sig_genes)}")
print(sig_genes)

# ============================================================
# 2. RECUPERER LES SYMBOLES DE GENES (depuis un des fichiers sources)
# ============================================================
ref_file = "data/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")

sig_genes = sig_genes.join(gene_symbols)
sig_genes.to_csv("outputs/interaction_genes_significant.csv")
print("\n--- Genes d'interaction avec symboles ---")
print(sig_genes[["Gene Symbol", "baseMean", "log2FoldChange", "pvalue", "padj"]])

# ============================================================
# 3. EXPRESSION (log2 CPM) DE CES GENES PAR ECHANTILLON
# ============================================================
counts = pd.read_csv("outputs/counts_filtered_for_deseq.csv", index_col=0, sep=";")  # samples x genes
library_sizes = counts.sum(axis=1)
cpm = counts.div(library_sizes, axis=0) * 1e6
log2_cpm = np.log2(cpm + 1)

meta_rows = []
for sample in counts.index:
    parts = sample.split("_")
    meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1]})
metadata = pd.DataFrame(meta_rows).set_index("sample").loc[counts.index]

# ============================================================
# 4. GRAPHIQUE : un panneau par gene d'interaction
# ============================================================
n_genes = len(sig_genes)
if n_genes > 0:
    fig, axes = plt.subplots(1, n_genes, figsize=(5 * n_genes, 5), squeeze=False)
    axes = axes[0]

    colors = {"Healthy": "tab:blue", "Mutant": "tab:red"}

    for ax, gene_id in zip(axes, sig_genes.index):
        symbol = sig_genes.loc[gene_id, "Gene Symbol"]
        if pd.isna(symbol):
            symbol = gene_id
        expr = log2_cpm[gene_id]
        plot_df = pd.DataFrame({"expr": expr}).join(metadata)

        positions = {"Female_Healthy": 0, "Female_Mutant": 1, "Male_Healthy": 2, "Male_Mutant": 3}
        for (sex, genotype), grp in plot_df.groupby(["sex", "genotype"]):
            pos = positions[f"{sex}_{genotype}"]
            jitter = np.random.uniform(-0.1, 0.1, size=len(grp))
            ax.scatter([pos + j for j in jitter], grp["expr"], color=colors[genotype], s=60)

        means = plot_df.groupby(["sex", "genotype"])["expr"].mean()
        ax.plot([0, 1], [means[("Female", "Healthy")], means[("Female", "Mutant")]],
                color="black", linestyle="--", alpha=0.6, label="Female")
        ax.plot([2, 3], [means[("Male", "Healthy")], means[("Male", "Mutant")]],
                color="gray", linestyle="--", alpha=0.6, label="Male")

        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(["F-Healthy", "F-Mutant", "M-Healthy", "M-Mutant"], rotation=45)
        ax.set_ylabel("log2(CPM+1)")
        ax.set_title(f"{symbol}\npadj={sig_genes.loc[gene_id, 'padj']:.4f}")

    plt.tight_layout()
    plt.savefig("outputs/interaction_genes_expression.pdf")
    plt.close()
    print("\nGraphique sauvegarde: outputs/interaction_genes_expression.png")
else:
    print("Aucun gene d'interaction significatif a tracer.")