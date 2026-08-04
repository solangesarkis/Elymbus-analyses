import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import pandas as pd
import gseapy as gp

os.makedirs("outputs_TG_final/kegg", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

ref_file = "data/Female Control NT TG.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ").str.upper()
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].replace("NAN", pd.NA)

print("Telechargement de la bibliotheque KEGG_2026 (souris)...")
kegg_sets = gp.get_library(name="KEGG_2026", organism="Mouse")

# Fond statistique = tous les genes reellement testes (pas le defaut Enrichr)
counts_for_deseq = read_csv_auto("outputs_TG_final/counts_filtered_for_deseq.csv")
tested_gene_ids = counts_for_deseq.columns.tolist()
background = gene_symbols.loc[gene_symbols.index.intersection(tested_gene_ids), "Gene Symbol"].dropna().unique().tolist()
print(f"Fond statistique (genes testes): {len(background)} genes")

# ============================================================
# PARTIE 1 : ORA (genes significatifs uniquement)
# ============================================================
ora_targets = {
    "Female_Mutant_BIMminusPF_DEGs": ("outputs_TG_final/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv", None),
    "Male_Healthy_BIMplusBAK_DEGs":  ("outputs_TG_final/DE_Male_Healthy_BIMplusBAK_vs_Untreated.csv", None),
    "Male_Healthy_BIMminusPF_DEGs":  ("outputs_TG_final/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv", None),
    "Interaction_sex_X_genotype":   ("outputs_TG_final/DE_interaction_sex_TMale_X_genotype_TMutant.csv", None),
    "Interaction_sex_X_BIMplusBAK": ("outputs_TG_final/DE_interaction_sex_TMale_X_treatment_TBIM+BAK.csv", "Fam220a"),
}

print("\n########## PARTIE 1: ORA KEGG (fond personnalise) ##########")
for label, (path, exclude_symbol) in ora_targets.items():
    if not os.path.exists(path):
        print(f"ATTENTION: fichier introuvable -> {path}")
        continue
    df = read_csv_auto(path)
    sig = df[df["padj"] < 0.05].copy()
    sig = sig.join(gene_symbols)
    if exclude_symbol:
        sig = sig[sig["Gene Symbol"] != exclude_symbol.upper()]
    gene_list = sig["Gene Symbol"].dropna().unique().tolist()
    print(f"\n=== {label}: {len(gene_list)} genes ===")

    if len(gene_list) < 3:
        print("Trop peu de genes pour un ORA fiable, ignore.")
        continue

    enr = gp.enrichr(gene_list=gene_list, gene_sets=kegg_sets, background=background,
                      outdir=None, cutoff=1.0)
    res = enr.results.sort_values("Adjusted P-value")
    res.to_csv(f"outputs_TG_final/kegg/ORA_{label}.csv", sep=";", index=False)

    sig_pathways = res[res["Adjusted P-value"] < 0.05]
    if len(sig_pathways) > 0:
        print(sig_pathways[["Term", "Overlap", "Adjusted P-value", "Genes"]].to_string())
    else:
        print("Aucune voie KEGG significative.")
        print("Top 5 (non significatif, reference):")
        print(res[["Term", "Overlap", "Adjusted P-value"]].head(5).to_string())

# ============================================================
# PARTIE 2 : GSEA preranked (tous les genes classes)
# ============================================================
def build_ranked_list(path):
    df = read_csv_auto(path)
    df = df.join(gene_symbols)
    df = df.dropna(subset=["Gene Symbol", "stat"])
    df["abs_stat"] = df["stat"].abs()
    df = df.sort_values("abs_stat", ascending=False).drop_duplicates("Gene Symbol")
    rnk = df[["Gene Symbol", "stat"]].sort_values("stat", ascending=False)
    rnk.columns = ["gene_name", "stat"]
    return rnk

gsea_targets = {
    "Female_Mutant_BIMminusPF": "outputs_TG_final/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv",
    "Male_Healthy_BIMplusBAK":  "outputs_TG_final/DE_Male_Healthy_BIMplusBAK_vs_Untreated.csv",
    "Male_Healthy_BIMminusPF":  "outputs_TG_final/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv",
    "Interaction_sex_X_genotype": "outputs_TG_final/DE_interaction_sex_TMale_X_genotype_TMutant.csv",
    "Interaction_sex_X_BIMplusBAK": "outputs_TG_final/DE_interaction_sex_TMale_X_treatment_TBIM+BAK.csv",
    "Interaction_sex_X_BIMminusPF": "outputs_TG_final/DE_interaction_sex_TMale_X_treatment_TBIM-PF.csv",
}

print("\n\n########## PARTIE 2: GSEA preranked KEGG ##########")
for label, path in gsea_targets.items():
    if not os.path.exists(path):
        print(f"ATTENTION: fichier introuvable -> {path}")
        continue
    rnk = build_ranked_list(path)
    print(f"\n=== {label}: {len(rnk)} genes classes ===")
    try:
        pre_res = gp.prerank(
            rnk=rnk, gene_sets=kegg_sets, permutation_num=1000,
            outdir=f"outputs_TG_final/kegg/GSEA_{label}", seed=42,
            min_size=5, max_size=1000,
        )
        res = pre_res.res2d
        res.to_csv(f"outputs_TG_final/kegg/GSEA_{label}_results.csv", sep=";")
        sig = res[res["FDR q-val"] < 0.05].sort_values("FDR q-val")
        print(f"{len(sig)} voies significatives (FDR<0.05) / {len(res)} testees")
        if len(sig) > 0:
            print(sig[["Term", "NES", "FDR q-val"]].head(10).to_string())
    except Exception as e:
        print(f"ERREUR: {e}")

print("\n=== TERMINE ===")
