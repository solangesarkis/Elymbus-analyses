import os
import pandas as pd
import gseapy as gp

os.makedirs("outputs_full", exist_ok=True)

# ============================================================
# 1. TELECHARGER LA BIBLIOTHEQUE GO BIOLOGICAL PROCESS (souris)
# ============================================================
names = gp.get_library_name(organism="Mouse")
go_bp_libs = [n for n in names if "go_biological_process" in n.lower()]
print("Bibliotheques GO_Biological_Process disponibles:", go_bp_libs)
go_lib_name = sorted(go_bp_libs)[-1]
print(f"Bibliotheque choisie: {go_lib_name}")

go_sets = gp.get_library(name=go_lib_name, organism="Mouse")
print(f"Nombre de voies GO BP: {len(go_sets)}")

# ============================================================
# 2. CHARGER ET NETTOYER LES 123 GENES SIGNIFICATIFS
# ============================================================
inter = pd.read_csv("outputs_full/interaction_triple_BIMPF_annotated.csv", sep=None, engine="python")
inter["Gene Symbol"] = inter["Gene Symbol"].astype(str).str.strip("'\" ")
gene_list_original = inter["Gene Symbol"].dropna().unique().tolist()
gene_list_upper = [g.upper() for g in gene_list_original]
print(f"\nNombre de genes uniques dans la liste: {len(gene_list_original)}")

all_go_genes = set()
for genes in go_sets.values():
    all_go_genes.update(genes)

overlap_original = len(set(gene_list_original) & all_go_genes)
overlap_upper = len(set(gene_list_upper) & all_go_genes)
print(f"Chevauchement casse originale: {overlap_original}/{len(gene_list_original)}")
print(f"Chevauchement casse majuscule: {overlap_upper}/{len(gene_list_original)}")

gene_list = gene_list_upper if overlap_upper > overlap_original else gene_list_original
print(f"Casse retenue: {'majuscules' if overlap_upper > overlap_original else 'originale'}")

# ============================================================
# 3. ANALYSE DE SUR-REPRESENTATION (ORA)
# ============================================================
enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets=go_sets,
    outdir="outputs_full/go_bp_enrichment_123genes",
    cutoff=1.0,
)

results = enr.results.sort_values("Adjusted P-value")
results.to_csv("outputs_full/go_bp_enrichment_123genes_results.csv", sep=";", index=False)

print("\n--- Voies significatives (Adjusted P-value < 0.05) ---")
pd.set_option("display.width", 250)
sig_results = results[results["Adjusted P-value"] < 0.05]
if len(sig_results) > 0:
    print(sig_results[["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]].to_string())
else:
    print("Aucune voie significative apres correction.")

print("\n--- Top 15 (toutes, meme non significatives, pour reference) ---")
print(results[["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]].head(15).to_string())

print(f"\n=== TERMINE. {len(sig_results)} voies significatives sur {len(results)} testees. ===")
print("Resultats complets: outputs_full/go_bp_enrichment_123genes_results.csv")
