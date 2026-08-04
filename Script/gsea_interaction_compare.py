import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import pandas as pd
import gseapy as gp

os.makedirs("outputs_full", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

# ============================================================
# 1. SYMBOLES DE GENES
# ============================================================
ref_file = "data/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ").str.upper()
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].replace("NAN", pd.NA)

def build_ranked_list(path):
    df = read_csv_auto(path)
    df = df.join(gene_symbols)
    df = df.dropna(subset=["Gene Symbol", "stat"])
    df["abs_stat"] = df["stat"].abs()
    df = df.sort_values("abs_stat", ascending=False).drop_duplicates("Gene Symbol")
    rnk = df[["Gene Symbol", "stat"]].sort_values("stat", ascending=False)
    rnk.columns = ["gene_name", "stat"]
    return rnk

# ============================================================
# 2. GSEA PRERANKED SUR L'INTERACTION TRIPLE POUR BIM+BAK
# ============================================================
path_bak = "outputs_full/DE_interaction_FIXED_sex_TMale_X_genotype_TMutant_X_treatment_TBIM+BAK.csv"
if not os.path.exists(path_bak):
    print(f"ATTENTION: fichier introuvable -> {path_bak}")
    print("Verifiez le nom exact avec: dir outputs_full\\DE_interaction_FIXED_sex*")
else:
    rnk_bak = build_ranked_list(path_bak)
    print(f"BIM+BAK: {len(rnk_bak)} genes classes")

    print("\nTelechargement de la bibliotheque GO_Biological_Process (souris)...")
    mouse_gene_sets = gp.get_library(name="GO_Biological_Process_2021", organism="Mouse")

    pre_res_bak = gp.prerank(
        rnk=rnk_bak,
        gene_sets=mouse_gene_sets,
        permutation_num=1000,
        outdir="outputs_full/gsea_interaction_triple_BIMplusBAK",
        seed=42,
        min_size=5,
        max_size=1000,
    )
    results_bak = pre_res_bak.res2d
    results_bak.to_csv("outputs_full/gsea_interaction_triple_BIMplusBAK_results.csv", sep=";")
    print("GSEA BIM+BAK termine.")

# ============================================================
# 3. COMPARER BIM+BAK vs BIM-PF (deja calcule precedemment)
# ============================================================
path_pf_results = "outputs_full/gsea_interaction_triple_BIMPF_results.csv"
path_bak_results = "outputs_full/gsea_interaction_triple_BIMplusBAK_results.csv"

if os.path.exists(path_pf_results) and os.path.exists(path_bak_results):
    res_pf = pd.read_csv(path_pf_results, sep=";")
    res_bak = pd.read_csv(path_bak_results, sep=";")

    pd.set_option("display.width", 250)

    print("\n\n========== TOP 15 VOIES - INTERACTION TRIPLE BIM+BAK ==========")
    top_bak = res_bak.sort_values("FDR q-val").head(15)
    print(top_bak[["Term", "NES", "NOM p-val", "FDR q-val"]].to_string())

    print("\n\n========== TOP 15 VOIES - INTERACTION TRIPLE BIM-PF ==========")
    top_pf = res_pf.sort_values("FDR q-val").head(15)
    print(top_pf[["Term", "NES", "NOM p-val", "FDR q-val"]].to_string())

    # Voies communes aux deux tops 15
    common = set(top_bak["Term"]) & set(top_pf["Term"])
    print(f"\n\nVoies communes entre les deux tops 15: {len(common)}")
    for t in common:
        print(" -", t)

    # Sauvegarder un tableau comparatif complet (toutes les voies significatives dans au moins un des deux)
    sig_bak = res_bak[res_bak["FDR q-val"] < 0.05][["Term", "NES", "FDR q-val"]].rename(
        columns={"NES": "NES_BIMplusBAK", "FDR q-val": "FDR_BIMplusBAK"})
    sig_pf = res_pf[res_pf["FDR q-val"] < 0.05][["Term", "NES", "FDR q-val"]].rename(
        columns={"NES": "NES_BIMminusPF", "FDR q-val": "FDR_BIMminusPF"})
    comparison = pd.merge(sig_bak, sig_pf, on="Term", how="outer")
    comparison.to_csv("outputs_full/gsea_interaction_comparison_BAK_vs_PF.csv", sep=";", index=False)
    print(f"\nTableau comparatif complet sauvegarde: outputs_full/gsea_interaction_comparison_BAK_vs_PF.csv")
else:
    print("\nATTENTION: un des deux fichiers de resultats GSEA est manquant.")
    print(f"BIM-PF existe: {os.path.exists(path_pf_results)}")
    print(f"BIM+BAK existe: {os.path.exists(path_bak_results)}")

print("\n=== TERMINE ===")
