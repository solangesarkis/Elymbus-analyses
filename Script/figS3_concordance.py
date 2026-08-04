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
from scipy.stats import pearsonr, spearmanr, binomtest

os.makedirs("outputs_full/figures", exist_ok=True)

def read_csv_auto(path, **kwargs):
    try:
        return pd.read_csv(path, sep=None, engine="python", encoding="utf-8", **kwargs)
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=None, engine="python", encoding="latin-1", **kwargs)

def load_lfc_padj(path):
    df = read_csv_auto(path).dropna(subset=["log2FoldChange", "padj"])
    return df[["log2FoldChange", "padj"]]

def concordance_panel(ax, path_x, path_y, xlabel, ylabel, panel_letter, title):
    dfx = load_lfc_padj(path_x).rename(columns={"log2FoldChange": "LFC_x", "padj": "padj_x"})
    dfy = load_lfc_padj(path_y).rename(columns={"log2FoldChange": "LFC_y", "padj": "padj_y"})
    merged = dfx.join(dfy, how="inner").dropna()

    sig = merged[(merged["padj_x"] < 0.05) | (merged["padj_y"] < 0.05)]
    n = len(sig)

    if n < 3:
        ax.text(0.5, 0.5, f"n={n}\n(too few for stats)", ha="center", va="center",
                transform=ax.transAxes, fontsize=7)
        ax.set_title(f"{panel_letter}  {title}", fontsize=7, fontweight="bold", loc="left")
        return

    x, y = sig["LFC_x"].values, sig["LFC_y"].values
    r, _ = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    slope = np.polyfit(x, y, 1)[0]
    same_dir = (np.sign(x) == np.sign(y))
    pct_same = 100 * same_dir.mean()

    ax.scatter(x, y, s=5, color="#8a9a9a", alpha=0.5, linewidths=0)
    xs_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs_line, np.poly1d(np.polyfit(x, y, 1))(xs_line), color="#b3123f", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)

    stats_text = (f"n={n}\nr={r:.2f}, rho={rho:.2f}\nslope={slope:.2f}\n{pct_same:.0f}% concordant")
    ax.text(0.03, 0.97, stats_text, transform=ax.transAxes, va="top", ha="left", fontsize=5.5,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", pad=0.2))

    ax.set_xlabel(xlabel, fontsize=6)
    ax.set_ylabel(ylabel, fontsize=6)
    ax.tick_params(axis="both", labelsize=5)
    ax.set_title(f"{panel_letter}  {title}", fontsize=7, fontweight="bold", loc="left")

D = "outputs_full"
panels = [
    ("A", f"{D}/DE_Female_Healthy_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv",
     "BIM+BAK LFC", "BIM-PF GEL LFC", "Product axis: Female Pitx2+/+"),
    ("B", f"{D}/DE_Female_Mutant_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv",
     "BIM+BAK LFC", "BIM-PF GEL LFC", "Product axis: Female Pitx2egl1/egl1"),
    ("C", f"{D}/DE_Male_Healthy_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv",
     "BIM+BAK LFC", "BIM-PF GEL LFC", "Product axis: Male Pitx2+/+"),
    ("D", f"{D}/DE_Male_Mutant_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv",
     "BIM+BAK LFC", "BIM-PF GEL LFC", "Product axis: Male Pitx2egl1/egl1"),

    ("E", f"{D}/DE_Female_Healthy_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Male_Healthy_BIMplusBAK_vs_Untreated.csv",
     "Female LFC", "Male LFC", "Sex axis: Pitx2+/+, BIM+BAK"),
    ("F", f"{D}/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv", f"{D}/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv",
     "Female LFC", "Male LFC", "Sex axis: Pitx2+/+, BIM-PF GEL"),
    ("G", f"{D}/DE_Female_Mutant_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Male_Mutant_BIMplusBAK_vs_Untreated.csv",
     "Female LFC", "Male LFC", "Sex axis: Pitx2egl1/egl1, BIM+BAK"),
    ("H", f"{D}/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv", f"{D}/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv",
     "Female LFC", "Male LFC", "Sex axis: Pitx2egl1/egl1, BIM-PF GEL"),

    ("I", f"{D}/DE_Female_Healthy_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Female_Mutant_BIMplusBAK_vs_Untreated.csv",
     "Pitx2+/+ LFC", "Pitx2egl1/egl1 LFC", "Genotype axis: Female, BIM+BAK"),
    ("J", f"{D}/DE_Female_Healthy_BIMminusPF_vs_Untreated.csv", f"{D}/DE_Female_Mutant_BIMminusPF_vs_Untreated.csv",
     "Pitx2+/+ LFC", "Pitx2egl1/egl1 LFC", "Genotype axis: Female, BIM-PF GEL"),
    ("K", f"{D}/DE_Male_Healthy_BIMplusBAK_vs_Untreated.csv", f"{D}/DE_Male_Mutant_BIMplusBAK_vs_Untreated.csv",
     "Pitx2+/+ LFC", "Pitx2egl1/egl1 LFC", "Genotype axis: Male, BIM+BAK"),
    ("L", f"{D}/DE_Male_Healthy_BIMminusPF_vs_Untreated.csv", f"{D}/DE_Male_Mutant_BIMminusPF_vs_Untreated.csv",
     "Pitx2+/+ LFC", "Pitx2egl1/egl1 LFC", "Genotype axis: Male, BIM-PF GEL"),
]

fig, axes = plt.subplots(3, 4, figsize=(8.27, 11.0))
axes_flat = axes.flatten()

for ax, (letter, path_x, path_y, xlabel, ylabel, title) in zip(axes_flat, panels):
    if not os.path.exists(path_x) or not os.path.exists(path_y):
        ax.text(0.5, 0.5, "File not found", ha="center", va="center", transform=ax.transAxes, fontsize=6)
        ax.set_title(f"{letter}  {title}", fontsize=7, fontweight="bold", loc="left")
        continue
    concordance_panel(ax, path_x, path_y, xlabel, ylabel, letter, title)

fig.suptitle("Gene-level concordance and directionality across product and host-context contrasts\n"
             "Gene selection: union of genes with padj<0.05 in either condition",
             fontsize=8, y=1.01)

plt.tight_layout()
plt.savefig("outputs_full/figures/FigureS3_concordance_all_axes_A4.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_full/figures/FigureS3_concordance_all_axes_A4.pdf")