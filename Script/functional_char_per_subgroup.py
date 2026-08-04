import os
import pandas as pd
import gseapy as gp

os.makedirs("outputs_full/functional_per_subgroup", exist_ok=True)

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

print("Telechargement des bibliotheques KEGG et GO BP (souris)...")
kegg_sets = gp.get_library(name="KEGG_2026", organism="Mouse")
go_sets = gp.get_library(name="GO_Biological_Process_2026", organism="Mouse")

def run_ora_and_annotate(gene_list_upper, gene_sets, df_full, label, db_name):
    """Lance l'ORA et retourne un tableau gene x voie avec direction up/down."""
    if len(gene_list_upper) < 3:
        print(f"  {db_name}: trop peu de genes ({len(gene_list_upper)}), ignore.")
        return None, None

    enr = gp.enrichr(gene_list=gene_list_upper, gene_sets=gene_sets, outdir=None, cutoff=1.0)
    res = enr.results.sort_values("Adjusted P-value")
    sig = res[res["Adjusted P-value"] < 0.05].copy()
    print(f"  {db_name}: {len(sig)} voies significatives / {len(res)} testees")

    if len(sig) == 0:
        return res, None

    # Construire le tableau gene x voie avec direction
    rows = []
    for _, row in sig.iterrows():
        term = row["Term"]
        padj_pathway = row["Adjusted P-value"]
        genes_in_term = row["Genes"].split(";")
        for g in genes_in_term:
            match = df_full[df_full["Gene Symbol"].str.upper() == g]
            if match.empty:
                continue
            lfc = match["log2FoldChange"].values[0]
            gene_padj = match["padj"].values[0]
            rows.append({
                "pathway": term, "pathway_padj": padj_pathway,
                "gene": g, "log2FC": lfc, "gene_padj": gene_padj,
                "direction": "up" if lfc > 0 else "down"
            })
    detail = pd.DataFrame(rows)
    return res, detail

all_summaries = []

for label, path in subgroups.items():
    print(f"\n========== {label} ==========")
    if not os.path.exists(path):
        print(f"ATTENTION: fichier introuvable -> {path}")
        continue

    df = read_csv_auto(path)
    df_full = df.join(gene_symbols)
    sig_df = df_full[df_full["padj"] < 0.05]
    gene_list = sig_df["Gene Symbol"].dropna().str.upper().unique().tolist()
    print(f"Genes significatifs (padj<0.05): {len(gene_list)}")

    # ---- KEGG ----
    kegg_res, kegg_detail = run_ora_and_annotate(gene_list, kegg_sets, df_full, label, "KEGG")
    if kegg_res is not None:
        kegg_res.to_csv(f"outputs_full/functional_per_subgroup/KEGG_all_{label}.csv", sep=";", index=False)
    if kegg_detail is not None:
        kegg_detail.to_csv(f"outputs_full/functional_per_subgroup/KEGG_genes_{label}.csv", sep=";", index=False)
        top_pathway = kegg_detail["pathway"].iloc[0] if not kegg_detail.empty else None
        if top_pathway:
            top_genes = kegg_detail[kegg_detail["pathway"] == top_pathway]
            n_up = (top_genes["direction"] == "up").sum()
            n_down = (top_genes["direction"] == "down").sum()
            print(f"  Top voie KEGG: {top_pathway} -> {n_up} up, {n_down} down")
            all_summaries.append({"subgroup": label, "db": "KEGG", "top_pathway": top_pathway,
                                   "n_up": n_up, "n_down": n_down,
                                   "padj": top_genes['pathway_padj'].iloc[0]})

    # ---- GO Biological Process ----
    go_res, go_detail = run_ora_and_annotate(gene_list, go_sets, df_full, label, "GO_BP")
    if go_res is not None:
        go_res.to_csv(f"outputs_full/functional_per_subgroup/GOBP_all_{label}.csv", sep=";", index=False)
    if go_detail is not None:
        go_detail.to_csv(f"outputs_full/functional_per_subgroup/GOBP_genes_{label}.csv", sep=";", index=False)
        top_pathway = go_detail["pathway"].iloc[0] if not go_detail.empty else None
        if top_pathway:
            top_genes = go_detail[go_detail["pathway"] == top_pathway]
            n_up = (top_genes["direction"] == "up").sum()
            n_down = (top_genes["direction"] == "down").sum()
            print(f"  Top voie GO BP: {top_pathway} -> {n_up} up, {n_down} down")
            all_summaries.append({"subgroup": label, "db": "GO_BP", "top_pathway": top_pathway,
                                   "n_up": n_up, "n_down": n_down,
                                   "padj": top_genes['pathway_padj'].iloc[0]})

summary_df = pd.DataFrame(all_summaries)
summary_df.to_csv("outputs_full/functional_per_subgroup/summary_top_pathways.csv", sep=";", index=False)
print("\n\n=== RESUME FINAL ===")
pd.set_option("display.width", 200)
print(summary_df.to_string(index=False))
print("\n=== TERMINE ===")
