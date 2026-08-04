import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

os.makedirs("outputs_full", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

# ============================================================
# 1. CHARGER LES VOIES GO SIGNIFICATIVES ET LEURS GENES
# ============================================================
go_results = pd.read_csv("outputs_full/go_bp_enrichment_123genes_results.csv", sep=";")
sig_terms = go_results[go_results["Adjusted P-value"] < 0.05].copy()
print(f"Voies significatives: {len(sig_terms)}")

gene_to_terms = {}
for _, row in sig_terms.iterrows():
    genes = row["Genes"].split(";")
    for g in genes:
        gene_to_terms.setdefault(g, []).append(row["Term"].split(" (GO:")[0])

all_genes_upper = sorted(gene_to_terms.keys())
print(f"Nombre de genes uniques impliques dans ces voies: {len(all_genes_upper)}")

# ============================================================
# 2. RETROUVER LES GENE ID (Ensembl) CORRESPONDANTS
# ============================================================
ref_file = "data/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ")
gene_symbols["Symbol_upper"] = gene_symbols["Gene Symbol"].str.upper()

symbol_to_id = gene_symbols.reset_index().drop_duplicates("Symbol_upper").set_index("Symbol_upper")["Gene ID"]

# ============================================================
# 3. RECUPERER LE LFC (BIM-PF vs Untreated) DANS LES 4 SOUS-GROUPES
# ============================================================
subgroup_files = {
    "Female_Healthy": "outputs_full/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv",
    "Female_Mutant":  "outputs_full/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv",
    "Male_Healthy":   "outputs_full/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv",
    "Male_Mutant":    "outputs_full/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv",
}

rows = []
for gene_upper in all_genes_upper:
    if gene_upper not in symbol_to_id.index:
        continue
    gene_id = symbol_to_id[gene_upper]
    row = {"Gene Symbol": gene_upper, "Pathways": " | ".join(sorted(set(gene_to_terms[gene_upper])))}
    for label, path in subgroup_files.items():
        if not os.path.exists(path):
            print(f"ATTENTION: fichier introuvable -> {path}")
            continue
        df = read_csv_auto(path)
        if gene_id in df.index:
            row[f"LFC_{label}"] = df.loc[gene_id, "log2FoldChange"]
            row[f"padj_{label}"] = df.loc[gene_id, "padj"]
    rows.append(row)

table = pd.DataFrame(rows).set_index("Gene Symbol")
table.to_csv("outputs_full/go_pathways_genes_across_conditions.csv")

lfc_cols = [c for c in table.columns if c.startswith("LFC_")]
print("\n--- LFC (BIM-PF vs Untreated) par sous-groupe, avec voie(s) associee(s) ---")
pd.set_option("display.width", 250)
print(table[["Pathways"] + lfc_cols].to_string())

# ============================================================
# 4. HEATMAP triee par voie biologique
# ============================================================
plot_df = table.dropna(subset=lfc_cols, how="all").copy()
plot_df = plot_df.sort_values("Pathways")

fig, ax = plt.subplots(figsize=(7, max(6, len(plot_df) * 0.25)))
mat = plot_df[lfc_cols].values.astype(float)
vmax = np.nanmax(np.abs(mat))
im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(4))
ax.set_xticklabels(["Female\nHealthy", "Female\nMutant", "Male\nHealthy", "Male\nMutant"])
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels([f"{g} ({p[:30]})" for g, p in zip(plot_df.index, plot_df["Pathways"])], fontsize=6)
ax.set_title("LFC (BIM-PF vs Untreated) pour les genes des voies GO significatives", fontsize=9)
cbar = plt.colorbar(im, ax=ax, shrink=0.5)
cbar.set_label("log2FoldChange")
plt.tight_layout()
plt.savefig("outputs_full/go_pathways_heatmap_by_subgroup.png", dpi=150)
plt.close()

print("\n=== TERMINE. Heatmap: outputs_full/go_pathways_heatmap_by_subgroup.png ===")
