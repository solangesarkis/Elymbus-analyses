import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

os.makedirs("outputs_full/figures", exist_ok=True)

path_pf = "outputs_full/gsea_interaction_triple_BIMPF_results.csv"
path_bak = "outputs_full/gsea_interaction_triple_BIMplusBAK_results.csv"

res_pf = pd.read_csv(path_pf, sep=None, engine="python")
res_bak = pd.read_csv(path_bak, sep=None, engine="python")

# Ne garder que les voies REELLEMENT significatives (FDR < 0.05)
res_pf_sig = res_pf[res_pf["FDR q-val"] < 0.05].copy()
res_bak_sig = res_bak[res_bak["FDR q-val"] < 0.05].copy()

res_pf_sig["treatment"] = "BIM-PF"
res_bak_sig["treatment"] = "BIM+BAK"

top_pf = res_pf_sig.sort_values("FDR q-val").head(10)
top_bak = res_bak_sig.sort_values("FDR q-val").head(10)
terms_of_interest = pd.concat([top_pf["Term"], top_bak["Term"]]).unique()

# Combiner, mais SEULEMENT les versions significatives (pas de point si non-significatif)
combined = pd.concat([res_pf_sig, res_bak_sig])
plot_df = combined[combined["Term"].isin(terms_of_interest)].copy()
plot_df["neglog10FDR"] = -np.log10(plot_df["FDR q-val"].clip(lower=1e-10))
plot_df["Term_short"] = plot_df["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True)

term_order = pd.concat([top_pf, top_bak])["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True).unique().tolist()
term_order = sorted(set(term_order), key=lambda t: plot_df[plot_df["Term_short"] == t]["NES"].mean())

fig, ax = plt.subplots(figsize=(8, max(6, len(term_order) * 0.35)))

for treat, xpos in [("BIM+BAK", 0), ("BIM-PF", 1)]:
    sub = plot_df[plot_df["treatment"] == treat]
    y = [term_order.index(t) for t in sub["Term_short"] if t in term_order]
    sub = sub[sub["Term_short"].isin(term_order)]
    sizes = sub["neglog10FDR"] * 40
    sc = ax.scatter([xpos] * len(sub), y, s=sizes, c=sub["NES"], cmap="RdBu_r",
                     vmin=-3, vmax=3, edgecolors="black", linewidths=0.5)

ax.set_yticks(range(len(term_order)))
ax.set_yticklabels(term_order, fontsize=8)
ax.set_xticks([0, 1])
ax.set_xticklabels(["BIM+BAK", "BIM-PF"])
ax.set_xlim(-0.5, 1.5)
ax.set_title("Top-10 interaction-term GSEA pathways per treatment\n(dot shown only if FDR < 0.05; size = -log10 FDR, color = NES)", fontsize=10)

cbar = plt.colorbar(sc, ax=ax, shrink=0.5)
cbar.set_label("NES")

plt.tight_layout()
plt.savefig("outputs_full/figures/fig4a_gsea_top20_fixed.pdf")
plt.close()
print("Sauvegarde: outputs_full/figures/fig4a_gsea_top20_fixed.pdf")
