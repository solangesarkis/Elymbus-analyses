import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import pandas as pd
import gseapy as gp

os.makedirs("outputs_full/gsea_per_subgroup", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

ref_file = "data/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ").str.upper()
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].replace("NAN", pd.NA)

subgroups = {
    "Female_Healthy": "outputs_full/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv",
    "Female_Mutant":  "outputs_full/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv",
    "Male_Healthy":   "outputs_full/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv",
    "Male_Mutant":    "outputs_full/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv",
}

def build_ranked_list(path):
    df = read_csv_auto(path)
    df = df.join(gene_symbols)
    df = df.dropna(subset=["Gene Symbol", "stat"])
    df["abs_stat"] = df["stat"].abs()
    df = df.sort_values("abs_stat", ascending=False).drop_duplicates("Gene Symbol")
    rnk = df[["Gene Symbol", "stat"]].sort_values("stat", ascending=False)
    rnk.columns = ["gene_name", "stat"]
    return rnk

print("Telechargement des bibliotheques (souris)...")
go_sets = gp.get_library(name="GO_Biological_Process_2026", organism="Mouse")
kegg_sets = gp.get_library(name="KEGG_2026", organism="Mouse")

all_top = []

for label, path in subgroups.items():
    print(f"\n========== {label} ==========")
    if not os.path.exists(path):
        print(f"ATTENTION: fichier introuvable -> {path}")
        continue

    rnk = build_ranked_list(path)
    print(f"Genes classes: {len(rnk)}")

    for db_name, gene_sets in [("GO_BP", go_sets), ("KEGG", kegg_sets)]:
        try:
            pre_res = gp.prerank(
                rnk=rnk,
                gene_sets=gene_sets,
                permutation_num=1000,
                outdir=f"outputs_full/gsea_per_subgroup/{label}_{db_name}",
                seed=42,
                min_size=5,
                max_size=1000,
            )
            res = pre_res.res2d
            res.to_csv(f"outputs_full/gsea_per_subgroup/{label}_{db_name}_results.csv", sep=";")
            sig = res[res["FDR q-val"] < 0.05].sort_values("FDR q-val")
            print(f"  {db_name}: {len(sig)} voies significatives (FDR<0.05) / {len(res)} testees")
            if len(sig) > 0:
                top = sig.iloc[0]
                print(f"    Top: {top['Term']} (NES={top['NES']:.2f}, FDR={top['FDR q-val']:.4g})")
                all_top.append({
                    "subgroup": label, "db": db_name, "n_significant": len(sig),
                    "top_term": top["Term"], "top_NES": top["NES"], "top_FDR": top["FDR q-val"]
                })
        except Exception as e:
            print(f"  {db_name}: ERREUR -> {e}")

summary_df = pd.DataFrame(all_top)
summary_df.to_csv("outputs_full/gsea_per_subgroup/summary_top_pathways.csv", sep=";", index=False)
print("\n\n=== RESUME FINAL ===")
pd.set_option("display.width", 200)
print(summary_df.to_string(index=False))
print("\n=== TERMINE ===")
