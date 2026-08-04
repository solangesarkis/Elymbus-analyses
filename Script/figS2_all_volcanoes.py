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

os.makedirs("outputs_full/figures", exist_ok=True)

def read_csv_auto(path, **kwargs):
    try:
        return pd.read_csv(path, sep=None, engine="python", encoding="utf-8", **kwargs)
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=None, engine="python", encoding="latin-1", **kwargs)

PADJ_THRESH = 0.05
LFC_THRESH = 1.0

comparisons = [
    ("A", "Female Pitx2+/+, BIM+BAK vs Untreated", "outputs_full/DE_Female_Healthy_BIMplusBAK_vs_Untreated.csv"),
    ("B", "Female Pitx2+/+, BIM-PF GEL vs Untreated", "outputs_full/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv"),
    ("C", "Male Pitx2+/+, BIM+BAK vs Untreated", "outputs_full/DE_Male_Healthy_BIMplusBAK_vs_Untreated.csv"),
    ("D", "Male Pitx2+/+, BIM-PF GEL vs Untreated", "outputs_full/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv"),
    ("E", "Female Pitx2egl1/egl1, BIM+BAK vs Untreated", "outputs_full/DE_Female_Mutant_BIMplusBAK_vs_Untreated.csv"),
    ("F", "Female Pitx2egl1/egl1, BIM-PF GEL vs Untreated", "outputs_full/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv"),
    ("G", "Male Pitx2egl1/egl1, BIM+BAK vs Untreated", "outputs_full/DE_Male_Mutant_BIMplusBAK_vs_Untreated.csv"),
    ("H", "Male Pitx2egl1/egl1, BIM-PF GEL vs Untreated", "outputs_full/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv"),
]

all_data = {}
global_max_lfc, global_max_neglog10padj = 0, 0

for letter, label, path in comparisons:
    if not os.path.exists(path):
        print(f"  ATTENTION [{letter}]: fichier introuvable -> {path}")
        continue
    df = read_csv_auto(path)
    df = df.dropna(subset=["padj", "log2FoldChange"])
    df["neglog10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    all_data[letter] = (label, df)
    global_max_lfc = max(global_max_lfc, df["log2FoldChange"].abs().max())
    global_max_neglog10padj = max(global_max_neglog10padj, df["neglog10padj"].max())

XLIM = global_max_lfc * 1.1
YLIM = global_max_neglog10padj * 1.1
print(f"\nLimites communes: X = +/-{XLIM:.2f}, Y = 0 a {YLIM:.2f}")

fig, axes = plt.subplots(4, 2, figsize=(8.0, 11.0), sharex=True, sharey=True)
axes_flat = axes.flatten()

POINT_SIZE_NOTSIG = 3
POINT_SIZE_SIG = 4.5

for ax, (letter, label, path) in zip(axes_flat, comparisons):
    if letter not in all_data:
        ax.text(0.5, 0.5, "Data not found", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"{letter}  {label}", fontsize=6.5, fontweight="bold", loc="left")
        continue

    _, df = all_data[letter]
    not_sig = df["padj"] >= PADJ_THRESH
    sig_up = (df["padj"] < PADJ_THRESH) & (df["log2FoldChange"] > 0)
    sig_down = (df["padj"] < PADJ_THRESH) & (df["log2FoldChange"] < 0)

    ax.scatter(df.loc[not_sig, "log2FoldChange"], df.loc[not_sig, "neglog10padj"],
               s=POINT_SIZE_NOTSIG, color="lightgray", alpha=0.5, linewidths=0)
    ax.scatter(df.loc[sig_up, "log2FoldChange"], df.loc[sig_up, "neglog10padj"],
               s=POINT_SIZE_SIG, color="#c0392b", alpha=0.8, linewidths=0)
    ax.scatter(df.loc[sig_down, "log2FoldChange"], df.loc[sig_down, "neglog10padj"],
               s=POINT_SIZE_SIG, color="#2166ac", alpha=0.8, linewidths=0)

    ax.axhline(-np.log10(PADJ_THRESH), color="black", linewidth=0.5, linestyle="--")
    ax.axvline(LFC_THRESH, color="black", linewidth=0.5, linestyle="--")
    ax.axvline(-LFC_THRESH, color="black", linewidth=0.5, linestyle="--")

    ax.set_xlim(-XLIM, XLIM)
    ax.set_ylim(0, YLIM)

    n_up, n_down = sig_up.sum(), sig_down.sum()
    ax.set_title(f"{letter}  {label}\n({n_up} up, {n_down} down)", fontsize=6.5, fontweight="bold", loc="left")
    ax.tick_params(axis="both", labelsize=5)

# Forcer l'affichage des nombres X sur TOUS les panneaux
for ax in axes_flat:
    ax.tick_params(axis="x", labelbottom=True, labelsize=5)

for ax in axes[-1, :]:
    ax.set_xlabel("log$_2$ fold change", fontsize=6)
for ax in axes[:, 0]:
    ax.set_ylabel("-log$_{10}$ adjusted p value", fontsize=6)

fig.suptitle("Complete volcano plots for all eight corneal treatment contrasts\n"
             "(common axes, thresholds: padj<0.05, |log2FC|>1)", fontsize=8, y=1.01)

plt.tight_layout()
plt.savefig("outputs_full/figures/FigureS2_all_volcanoes_A4.pdf", bbox_inches="tight")
plt.close()
print("\nSauvegarde: outputs_full/figures/FigureS2_all_volcanoes_A4.pdf")