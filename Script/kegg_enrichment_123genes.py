import os
import pandas as pd
import gseapy as gp

os.makedirs("outputs_full", exist_ok=True)

# ============================================================
# 1. TROUVER AUTOMATIQUEMENT LA BONNE BIBLIOTHEQUE KEGG (souris)
# ============================================================
names = gp.get_library_name(organism="Mouse")
kegg_libs = [n for n in names if "kegg" in n.lower()]
print("Bibliotheques KEGG disponibles:", kegg_libs)
kegg_lib_name = sorted(kegg_libs)[-1]  # prend la plus recente (ordre alphabetique/annee)
print(f"Bibliotheque choisie: {kegg_lib_name}")

kegg_sets = gp.get_library(name=kegg_lib_name, organism="Mouse")
print(f"Nombre de voies KEGG: {len(kegg_sets)}")

# Verifier la casse utilisee dans la bibliotheque
sample_genes = list(kegg_sets.values())[0][:5]
print(f"Exemple de genes dans la bibliotheque (verifier la casse): {sample_genes}")

# ============================================================
# 2. CHARGER ET NETTOYER LES 123 GENES SIGNIFICATIFS
# ============================================================
inter = pd.read_csv("outputs_full/interaction_triple_BIMPF_annotated.csv", sep=None, engine="python")
inter["Gene Symbol"] = inter["Gene Symbol"].astype(str).str.strip("'\" ")
gene_list_original = inter["Gene Symbol"].dropna().unique().tolist()
gene_list_upper = [g.upper() for g in gene_list_original]
print(f"\nNombre de genes uniques dans la liste: {len(gene_list_original)}")

# Tester quelle casse donne le plus grand chevauchement avec la bibliotheque
all_kegg_genes = set()
for genes in kegg_sets.values():
    all_kegg_genes.update(genes)

overlap_original = len(set(gene_list_original) & all_kegg_genes)
overlap_upper = len(set(gene_list_upper) & all_kegg_genes)
print(f"Chevauchement casse originale: {overlap_original}/{len(gene_list_original)}")
print(f"Chevauchement casse majuscule: {overlap_upper}/{len(gene_list_original)}")

gene_list = gene_list_upper if overlap_upper > overlap_original else gene_list_original
print(f"Casse retenue: {'majuscules' if overlap_upper > overlap_original else 'originale'}")

# ============================================================
# 3. ANALYSE DE SUR-REPRESENTATION (ORA) KEGG
# ============================================================
enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets=kegg_sets,
    outdir="outputs_full/kegg_enrichment_123genes",
    cutoff=1.0,  # on garde tout, on filtrera nous-memes
)

results = enr.results.sort_values("Adjusted P-value")
results.to_csv("outputs_full/kegg_enrichment_123genes_results.csv", sep=";", index=False)

print("\n--- Top 15 voies KEGG (triees par p-value ajustee) ---")
pd.set_option("display.width", 200)
print(results[["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]].head(15).to_string())

print("\n=== TERMINE. Resultats: outputs_full/kegg_enrichment_123genes_results.csv ===")
