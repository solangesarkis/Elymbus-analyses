import os
import pandas as pd

os.makedirs("outputs_TG_full", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

# ============================================================
# 1. GENES SUSPECTS ET LEURS ID
# ============================================================
suspect_symbols = ["Krt12", "Prl", "Pomc", "Gh", "Cga", "Lrrc30", "Hjv"]

ref_file = "data/Female Control NT TG.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ")
suspect_ids = gene_symbols[gene_symbols["Gene Symbol"].isin(suspect_symbols)]
print("Gene IDs trouves:")
print(suspect_ids)

# ============================================================
# 2. RECUPERER LES COMPTAGES BRUTS (non normalises) POUR CES GENES
# ============================================================
counts_rounded = read_csv_auto("outputs_TG_full/counts_rounded.csv")  # genes x samples

meta_rows = []
for sample in counts_rounded.columns:
    parts = sample.split("_")
    meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1], "treatment": parts[2]})
metadata = pd.DataFrame(meta_rows).set_index("sample")

print("\n--- Comptages bruts (arrondis) pour les genes suspects, par echantillon ---")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", None)

for gene_id, symbol in zip(suspect_ids.index, suspect_ids["Gene Symbol"]):
    if gene_id not in counts_rounded.index:
        print(f"\n{symbol} ({gene_id}): absent apres filtrage (deja tres faible partout)")
        continue
    expr = counts_rounded.loc[gene_id]
    expr_df = pd.DataFrame({"count": expr}).join(metadata)
    expr_df = expr_df.sort_values("count", ascending=False)
    n_nonzero = (expr_df["count"] > 0).sum()
    print(f"\n=== {symbol} ({gene_id}) — {n_nonzero}/60 echantillons non-nuls ===")
    print(expr_df.head(10).to_string())

print("\n=== TERMINE ===")
