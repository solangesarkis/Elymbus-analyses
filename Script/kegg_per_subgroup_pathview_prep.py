import os
import pandas as pd
import requests
import gseapy as gp

os.makedirs("outputs_full/pathview", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

ref_file = "data/Cornea Read counts Female Healthy Untreated 1.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ")

subgroups = {
    "Female_Healthy": "outputs_full/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv",
    "Female_Mutant":  "outputs_full/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv",
    "Male_Healthy":   "outputs_full/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv",
    "Male_Mutant":    "outputs_full/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv",
}

# ============================================================
# 1. RECUPERER LA LISTE DES VOIES KEGG SOURIS (nom -> ID mmuXXXXX)
#    via l'API KEGG officielle
# ============================================================
print("Telechargement de la liste des voies KEGG souris (mmu)...")
resp = requests.get("https://rest.kegg.org/list/pathway/mmu")
kegg_id_map = {}
for line in resp.text.strip().split("\n"):
    kegg_id, name = line.split("\t")
    kegg_id = kegg_id.replace("path:", "")
    name_clean = name.split(" - Mus musculus")[0].strip()
    kegg_id_map[name_clean.upper()] = kegg_id
print(f"{len(kegg_id_map)} voies KEGG souris chargees.")

# ============================================================
# 2. KEGG (bibliotheque mouse) pour l'ORA
# ============================================================
kegg_sets = gp.get_library(name="KEGG_2026", organism="Mouse")

summary_rows = []

for label, path in subgroups.items():
    if not os.path.exists(path):
        print(f"ATTENTION: fichier introuvable -> {path}")
        continue

    df = read_csv_auto(path)
    df_full = df.join(gene_symbols)  # pour l'export complet (up ET down)

    # --- ORA sur les genes significatifs ---
    sig = df_full[df_full["padj"] < 0.05]
    gene_list = sig["Gene Symbol"].dropna().str.upper().unique().tolist()
    print(f"\n=== {label}: {len(gene_list)} genes significatifs (padj<0.05) ===")

    if len(gene_list) < 3:
        print("Trop peu de genes pour un ORA fiable, ignore.")
        continue

    enr = gp.enrichr(gene_list=gene_list, gene_sets=kegg_sets, outdir=None, cutoff=1.0)
    res = enr.results.sort_values("Adjusted P-value")
    top = res.iloc[0]
    top_term = top["Term"]
    top_padj = top["Adjusted P-value"]
    print(f"Top voie KEGG: {top_term} (padj={top_padj:.4g})")

    # Chercher l'ID KEGG correspondant
    kegg_id = kegg_id_map.get(top_term.upper())
    if kegg_id is None:
        # tentative de correspondance partielle
        matches = [v for k, v in kegg_id_map.items() if top_term.upper() in k or k in top_term.upper()]
        kegg_id = matches[0] if matches else None
    print(f"ID KEGG trouve: {kegg_id}")

    summary_rows.append({
        "subgroup": label, "top_KEGG_term": top_term,
        "adjusted_pvalue": top_padj, "kegg_id": kegg_id,
        "n_sig_genes": len(gene_list)
    })

    # --- Export TOUS les genes testes (up et down, pas que les significatifs)
    #     avec leur log2FC, pour colorer la pathview map ---
    export = df_full.dropna(subset=["Gene Symbol", "log2FoldChange"])
    export = export.drop_duplicates("Gene Symbol")[["Gene Symbol", "log2FoldChange"]]
    export.columns = ["symbol", "log2FC"]
    out_path = f"outputs_full/pathview/{label}_log2FC_for_pathview.csv"
    export.to_csv(out_path, index=False)
    print(f"Export log2FC complet: {out_path} ({len(export)} genes)")

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("outputs_full/pathview/top_kegg_per_subgroup.csv", index=False)
print("\n=== RESUME ===")
print(summary_df.to_string(index=False))
print("\n=== TERMINE. Verifiez outputs_full/pathview/top_kegg_per_subgroup.csv avant de lancer R ===")
