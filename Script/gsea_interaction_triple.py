import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import pandas as pd
import gseapy as gp

os.makedirs("outputs_full", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

# ============================================================
# 1. RECUPERER ET NETTOYER LES SYMBOLES DE GENES
# ============================================================
ref_file = "data/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ").str.upper()
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].replace("NAN", pd.NA)

# ============================================================
# 2. CONSTRUIRE LA LISTE CLASSEE (tous les genes testes, par stat)
# ============================================================
path = "outputs_full/DE_interaction_FIXED_sex_TMale_X_genotype_TMutant_X_treatment_TBIM-PF.csv"
df = read_csv_auto(path)
print(f"Genes testes dans l'interaction: {len(df)}")

df = df.join(gene_symbols)
df = df.dropna(subset=["Gene Symbol", "stat"])
df["abs_stat"] = df["stat"].abs()
df = df.sort_values("abs_stat", ascending=False).drop_duplicates("Gene Symbol")
rnk = df[["Gene Symbol", "stat"]].sort_values("stat", ascending=False)
rnk.columns = ["gene_name", "stat"]
rnk.to_csv("outputs_full/interaction_triple_BIMPF.rnk", sep="\t", index=False, header=False)
print(f"Genes classes pour GSEA: {len(rnk)}")

# ============================================================
# 3. TELECHARGER LA BIBLIOTHEQUE SOURIS ET LANCER GSEA
# ============================================================
print("\nTelechargement de la bibliotheque GO_Biological_Process_2021 (souris)...")
mouse_gene_sets = gp.get_library(name="GO_Biological_Process_2021", organism="Mouse")

pre_res = gp.prerank(
    rnk=rnk,
    gene_sets=mouse_gene_sets,
    permutation_num=1000,
    outdir="outputs_full/gsea_interaction_triple_BIMPF",
    seed=42,
    min_size=5,
    max_size=1000,
)
results = pre_res.res2d
results.to_csv("outputs_full/gsea_interaction_triple_BIMPF_results.csv", sep=";")

fdr_col = [c for c in results.columns if "fdr" in c.lower() or "q-val" in c.lower()]
results_sorted = results.sort_values(fdr_col[0]) if fdr_col else results
print("\nTop 15 voies - interaction triple (BIM-PF):")
print(results_sorted[["Term", "NES", "NOM p-val", "FDR q-val"]].head(15).to_string())

print("\n=== TERMINE ===")
