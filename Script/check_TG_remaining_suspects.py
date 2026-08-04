import os
import pandas as pd

os.makedirs("outputs_TG_full_clean", exist_ok=True)

def read_csv_auto(path):
    return pd.read_csv(path, index_col=0, sep=None, engine="python")

# ============================================================
# 1. GENES SUSPECTS
# ============================================================
suspect_symbols = ["Krt12", "Hjv", "Acta1", "Ckm", "Mb", "Tnnt3", "Myl1",
                    "Myh1", "Tnni2", "Ryr1", "Ttn", "Tnnc2", "Mylpf"]

ref_file = "data/Female Control NT TG.xlsx"
gene_symbols = pd.read_excel(ref_file, usecols=["Gene ID", "Gene Symbol"]).set_index("Gene ID")
gene_symbols["Gene Symbol"] = gene_symbols["Gene Symbol"].astype(str).str.strip("'\" ")
suspect_ids = gene_symbols[gene_symbols["Gene Symbol"].isin(suspect_symbols)]

# ============================================================
# 2. CHARGER LES COMPTAGES (donnees nettoyees, 58 echantillons)
# ============================================================
counts_for_deseq = read_csv_auto("outputs_TG_full_clean/counts_filtered_for_deseq.csv")  # samples x genes
counts_rounded = counts_for_deseq.T  # transposer en genes x samples

meta_rows = []
for sample in counts_rounded.columns:
    parts = sample.split("_")
    meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1], "treatment": parts[2]})
metadata = pd.DataFrame(meta_rows).set_index("sample")

pd.set_option("display.width", 250)
print("--- Comptages (donnees nettoyees) pour les genes suspects restants ---")

for gene_id, symbol in zip(suspect_ids.index, suspect_ids["Gene Symbol"]):
    if gene_id not in counts_rounded.index:
        print(f"\n{symbol} ({gene_id}): absent")
        continue
    expr = counts_rounded.loc[gene_id]
    expr_df = pd.DataFrame({"count": expr}).join(metadata)
    expr_df = expr_df.sort_values("count", ascending=False)
    n_nonzero = (expr_df["count"] > 0).sum()
    print(f"\n=== {symbol} ({gene_id}) — {n_nonzero}/58 echantillons non-nuls ===")
    print(expr_df.head(10).to_string())

print("\n=== TERMINE ===")
