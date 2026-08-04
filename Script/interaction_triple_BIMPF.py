import os
import pandas as pd

os.makedirs("outputs_full", exist_ok=True)

# ============================================================
# 1. CHARGER LES 123 GENES D'INTERACTION TRIPLE (BIM-PF)
# ============================================================
inter_path = "outputs_full/DE_interaction_FIXED_sex_TMale_X_genotype_TMutant_X_treatment_TBIM-PF.csv"
inter = pd.read_csv(inter_path, index_col=0, sep=None, engine="python")
sig = inter[inter["padj"] < 0.05].sort_values("padj")
print(f"Nombre de genes d'interaction triple (BIM-PF): {len(sig)}")

# ============================================================
# 2. AJOUTER LES SYMBOLES DE GENES
# ============================================================
ref_file = "data/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ")
sig = sig.join(gene_symbols)

# ============================================================
# 3. RECUPERER LE log2FC (traitement vs Untreated) DANS CHAQUE SOUS-GROUPE
#    (fichiers deja generes par la Partie A)
# ============================================================
subgroup_files = {
    "Female_Healthy": "outputs_full/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv",
    "Female_Mutant":  "outputs_full/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv",
    "Male_Healthy":   "outputs_full/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv",
    "Male_Mutant":    "outputs_full/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv",
}

for label, path in subgroup_files.items():
    df = pd.read_csv(path, index_col=0, sep=None, engine="python")
    sig[f"LFC_{label}"] = df.reindex(sig.index)["log2FoldChange"]
    sig[f"padj_{label}"] = df.reindex(sig.index)["padj"]

# ============================================================
# 4. SAUVEGARDER LE TABLEAU COMPLET
# ============================================================
out_cols = ["Gene Symbol", "log2FoldChange", "padj"] + \
           [c for c in sig.columns if c.startswith("LFC_") or c.startswith("padj_")]
final = sig[out_cols]
final.to_csv("outputs_full/interaction_triple_BIMPF_annotated.csv")

print("\n--- Top 20 genes (avec LFC par sous-groupe pour verifier la coherence de direction) ---")
pd.set_option("display.width", 200)
print(final.head(20).to_string())

# Compter combien de genes ont des directions incoherentes entre sous-groupes
lfc_cols = [c for c in final.columns if c.startswith("LFC_")]
signs = final[lfc_cols].apply(lambda row: set(x > 0 for x in row.dropna()), axis=1)
inconsistent = signs.apply(lambda s: len(s) > 1)
print(f"\nGenes avec directions incoherentes entre les 4 sous-groupes: {inconsistent.sum()} / {len(final)}")

print("\n=== TERMINE. Tableau complet: outputs_full/interaction_triple_BIMPF_annotated.csv ===")
