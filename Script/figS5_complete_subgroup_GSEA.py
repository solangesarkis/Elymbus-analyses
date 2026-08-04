import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

os.makedirs("outputs_full/figures", exist_ok=True)

def read_csv_auto(path, sep=None):
    try:
        return pd.read_csv(path, sep=sep, engine="python", encoding="utf-8")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=sep, engine="python", encoding="latin-1")

NES_VMIN, NES_VMAX = -3, 3
CMAP = "RdBu_r"
FDR_THRESH = 0.05
JACCARD_THRESHOLD = 0.5
MIN_READABLE_FONTSIZE = 6.0       # plancher de lisibilite absolu
MAX_TERMS_ONE_PAGE = 55           # au-dela, meme a 6pt, ca ne tient plus sur une A4

subgroups_bak = {s: f"outputs_full/gsea_per_subgroup_BIMBAK/gsea_per-subgroup_BIMBAK/{s}_GO_BP_results.csv"
                  for s in ["Female_Healthy", "Male_Healthy", "Female_Mutant", "Male_Mutant"]}
subgroups_pf = {s: f"outputs_full/gsea_per_subgroup/{s}_GO_BP_results.csv"
                 for s in ["Female_Healthy", "Male_Healthy", "Female_Mutant", "Male_Mutant"]}

col_labels_display = ["Female\nPitx2+/+", "Male\nPitx2+/+", "Female\nPitx2egl1/egl1", "Male\nPitx2egl1/egl1"]

def load_all_subgroups(subgroup_files):
    all_data = {}
    for sub, path in subgroup_files.items():
        if not os.path.exists(path):
            print(f"  ATTENTION: introuvable -> {path}")
            continue
        df = read_csv_auto(path)
        df["Term_short"] = df["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True)
        all_data[sub] = df.set_index("Term_short")
    return all_data

data_bak = load_all_subgroups(subgroups_bak)
data_pf = load_all_subgroups(subgroups_pf)
data_combined = {**data_bak, **data_pf}

def mean_nes_for_ordering(term):
    values = []
    for data_dict in [data_bak, data_pf]:
        for sub, df in data_dict.items():
            if term in df.index:
                row = df.loc[term]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                if row["FDR q-val"] < FDR_THRESH:
                    values.append(row["NES"])
    return np.mean(values) if values else 0

def get_lead_genes(term, data_dict):
    for sub, df in data_dict.items():
        if term in df.index:
            row = df.loc[term]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            if row["FDR q-val"] < FDR_THRESH and "Lead_genes" in row:
                return set(str(row["Lead_genes"]).split(";"))
    return set()

all_sig_terms = set()
for data_dict in [data_bak, data_pf]:
    for sub, df in data_dict.items():
        sig_terms = df[df["FDR q-val"] < FDR_THRESH].index
        all_sig_terms.update(sig_terms)

term_best_fdr = {}
for term in all_sig_terms:
    best = 1.0
    for data_dict in [data_bak, data_pf]:
        for sub, df in data_dict.items():
            if term in df.index:
                row = df.loc[term]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                best = min(best, row["FDR q-val"])
    term_best_fdr[term] = best

sorted_terms = sorted(all_sig_terms, key=lambda t: term_best_fdr[t])

kept_terms = []
kept_gene_sets = []
for term in sorted_terms:
    genes = get_lead_genes(term, data_combined)
    is_redundant = False
    for kept_genes in kept_gene_sets:
        union = genes | kept_genes
        inter = genes & kept_genes
        jaccard = len(inter) / len(union) if len(union) > 0 else 0
        if jaccard > JACCARD_THRESHOLD:
            is_redundant = True
            break
    if not is_redundant:
        kept_terms.append(term)
        kept_gene_sets.append(genes)

print(f"Voies avant dedoublonnage: {len(all_sig_terms)}")
print(f"Voies apres dedoublonnage: {len(kept_terms)}")

# ============================================================
# GARANTIR UNE SEULE PAGE : si trop de voies meme apres dedoublonnage,
# garder seulement les MAX_TERMS_ONE_PAGE au signal le plus fort (meilleur FDR)
# ============================================================
n_dropped = 0
if len(kept_terms) > MAX_TERMS_ONE_PAGE:
    kept_terms_sorted_by_signal = sorted(kept_terms, key=lambda t: term_best_fdr[t])
    n_dropped = len(kept_terms) - MAX_TERMS_ONE_PAGE
    kept_terms = kept_terms_sorted_by_signal[:MAX_TERMS_ONE_PAGE]
    print(f"ATTENTION: {n_dropped} voies supplementaires retirees pour tenir sur une page "
          f"(gardees: les {MAX_TERMS_ONE_PAGE} au meilleur FDR)")

term_order = sorted(kept_terms, key=mean_nes_for_ordering)
n_terms = len(term_order)

# Police calculee pour remplir exactement la page A4 disponible (hauteur ~9.5 pouces utiles)
label_fontsize = max(MIN_READABLE_FONTSIZE, min(9, 500 / n_terms))
print(f"Voies affichees: {n_terms}, police des labels: {label_fontsize:.1f}pt")

def plot_dotplot(ax, data_dict, subgroup_order, term_order, title):
    sc = None
    for j, sub in enumerate(subgroup_order):
        if sub not in data_dict:
            continue
        df = data_dict[sub]
        for i, term in enumerate(term_order):
            if term not in df.index:
                continue
            row = df.loc[term]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            if row["FDR q-val"] >= FDR_THRESH:
                continue
            size = max(-np.log10(max(row["FDR q-val"], 1e-10)) * 18, 8)
            sc = ax.scatter(j, i, s=size, c=row["NES"], cmap=CMAP, vmin=NES_VMIN, vmax=NES_VMAX,
                             edgecolors="black", linewidths=0.3)
    ax.set_xticks(range(len(subgroup_order)))
    ax.set_xticklabels(col_labels_display, fontsize=7)
    ax.set_yticks(range(len(term_order)))
    ax.set_yticklabels(term_order, fontsize=label_fontsize)
    ax.set_ylim(-0.7, len(term_order) - 0.3)
    ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
    return sc

subgroup_order = ["Female_Healthy", "Male_Healthy", "Female_Mutant", "Male_Mutant"]

fig = plt.figure(figsize=(8.27, 11.0))
gs = gridspec.GridSpec(1, 2, wspace=0.9, left=0.35, right=0.86, top=0.93, bottom=0.05)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])

plot_dotplot(ax1, data_bak, subgroup_order, term_order, "BIM+BAK vs Untreated")
sc_final = plot_dotplot(ax2, data_pf, subgroup_order, term_order, "BIM-PF GEL vs Untreated")
ax2.set_yticks([])

if sc_final is not None:
    cbar_ax = fig.add_axes([0.89, 0.55, 0.02, 0.32])
    cbar = fig.colorbar(sc_final, cax=cbar_ax)
    cbar.set_label("NES", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

legend_fdr_values = [0.05, 0.01, 0.001]
legend_sizes = [max(-np.log10(v) * 18, 8) for v in legend_fdr_values]
legend_ax = fig.add_axes([0.89, 0.15, 0.09, 0.22])
legend_ax.axis("off")
for i, (fdr_val, size) in enumerate(zip(legend_fdr_values, legend_sizes)):
    legend_ax.scatter([0.2], [1 - i * 0.35], s=size, c="gray", edgecolors="black", linewidths=0.3)
    legend_ax.text(0.5, 1 - i * 0.35, f"FDR={fdr_val}", fontsize=6, va="center")
legend_ax.set_xlim(0, 1.5)
legend_ax.set_ylim(-0.1, 1.2)
legend_ax.set_title("Significance", fontsize=6.5, loc="left")

subtitle = "(redundant pathways collapsed by gene-set overlap; FDR<0.05"
if n_dropped > 0:
    subtitle += f"; top {n_terms} strongest shown, {n_dropped} additional significant pathways omitted for space)"
else:
    subtitle += ")"

fig.suptitle(f"Complete direct subgroup pathway enrichment\n{subtitle}", fontsize=8, y=0.99)

plt.savefig("outputs_full/figures/FigureS5_complete_subgroup_GSEA_onepage.pdf", bbox_inches="tight")
plt.close()
print("\nSauvegarde: outputs_full/figures/FigureS5_complete_subgroup_GSEA_onepage.pdf")