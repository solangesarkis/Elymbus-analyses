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
from sklearn.decomposition import PCA

os.makedirs("outputs_TG_final/figures", exist_ok=True)

def read_csv_auto(path, **kwargs):
    try:
        return pd.read_csv(path, sep=None, engine="python", encoding="utf-8", **kwargs)
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=None, engine="python", encoding="latin-1", **kwargs)

def safe_load(path, label, **kwargs):
    if not os.path.exists(path):
        print(f"  ATTENTION [{label}]: fichier introuvable -> {path}")
        return None
    return read_csv_auto(path, **kwargs)

GENO_DISPLAY = {"Healthy": "Pitx2+/+", "Mutant": "Pitx2egl1/egl1"}
GENO_DISPLAY_MATH = {"Healthy": "Pitx2$^{+/+}$", "Mutant": "Pitx2$^{egl1/egl1}$"}
PADJ_THRESH = 0.05
LFC_THRESH = 1.0

fig, axes = plt.subplots(4, 2, figsize=(8.0, 11.0))
axA, axB, axC, axD, axE1, axE2, axF, axG = axes.flatten()

# ============================================================
# PANEL A
# ============================================================
scores = safe_load("outputs_TG_final/contamination_scores_with_flags.csv", "A - contamination")
if scores is not None:
    for i, col in enumerate(["muscle_score", "pituitary_score", "cornea_score"]):
        threshold = scores[col].mean() + 2 * scores[col].std()
        axA.scatter([i] * len(scores), scores[col], c=scores["flag_any"].map({True: "#d62728", False: "#2166ac"}),
                    s=14, edgecolors="black", linewidths=0.3, alpha=0.7)
        axA.hlines(threshold, i - 0.3, i + 0.3, color="red", linestyle="--", linewidth=1)
    axA.set_xticks([0, 1, 2])
    axA.set_xticklabels(["Muscle", "Pituitary", "Corneal"], fontsize=6)
    axA.set_ylabel("Composite score", fontsize=6)
    axA.tick_params(axis="y", labelsize=5)
    axA.set_title("A", fontsize=9, fontweight="bold", loc="left")

# ============================================================
# PANEL B + C
# ============================================================
counts_filtered = safe_load("outputs_TG_final/counts_filtered_for_deseq.csv", "B/C - counts", index_col=0)
metadata_tg = None
if counts_filtered is not None:
    meta_rows = []
    for sample in counts_filtered.index:
        parts = sample.split("_")
        meta_rows.append({"sample": sample, "sex": parts[0], "genotype": parts[1], "treatment": parts[2]})
    metadata_tg = pd.DataFrame(meta_rows).set_index("sample")
    retained_counts = metadata_tg.groupby(["sex", "genotype", "treatment"]).size().reset_index(name="n")
    labels_b = [f"{r.sex[0]} {GENO_DISPLAY_MATH[r.genotype]}\n{r.treatment}" for r in retained_counts.itertuples()]
    colors_b = ["#c0392b" if n < 3 else "#2166ac" for n in retained_counts["n"]]
    axB.bar(labels_b, retained_counts["n"], color=colors_b, edgecolor="black")
    axB.axhline(3, color="gray", linestyle=":", linewidth=0.8)
    axB.set_ylabel("Retained animals (n)", fontsize=6)
    axB.tick_params(axis="x", labelsize=4.2, rotation=90)
    axB.tick_params(axis="y", labelsize=5)
    axB.set_title("B", fontsize=9, fontweight="bold", loc="left")

    library_sizes = counts_filtered.sum(axis=1)
    cpm = counts_filtered.div(library_sizes, axis=0) * 1e6
    log2cpm = np.log2(cpm + 1)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(log2cpm)
    var_exp = pca.explained_variance_ratio_ * 100
    markers = {"Female": "o", "Male": "^"}
    colors_pca = {"Healthy": "#2166ac", "Mutant": "#d62728"}
    for (sex, genotype), grp in metadata_tg.groupby(["sex", "genotype"]):
        idx = [log2cpm.index.get_loc(s) for s in grp.index]
        axC.scatter(pcs[idx, 0], pcs[idx, 1], c=colors_pca[genotype], marker=markers[sex],
                    s=18, label=f"{sex} {GENO_DISPLAY[genotype]}", edgecolors="black", linewidths=0.3)
    axC.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)", fontsize=6)
    axC.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)", fontsize=6)
    axC.tick_params(axis="both", labelsize=5)
    axC.legend(fontsize=4.5)
    axC.set_title("C", fontsize=9, fontweight="bold", loc="left")

# ============================================================
# PANEL D
# ============================================================
de_summary = safe_load("outputs_TG_final/Table_TG_treatment_DEG_summary.csv", "D - DE summary")
if de_summary is not None:
    label_cols = [c for c in de_summary.columns if de_summary[c].dtype == object]
    n_col = [c for c in de_summary.columns if "n_DEG" in c or "n_sig" in c.lower()]
    n_col = n_col[0] if n_col else de_summary.select_dtypes(include="number").columns[0]
    labels_d = []
    for _, r in de_summary.iterrows():
        parts = []
        for c in label_cols:
            val = str(r[c])
            val = val.replace("Mutant", "Pitx2egl1/egl1").replace("Healthy", "Pitx2+/+")
            parts.append(val)
        labels_d.append(" ".join(parts))
    axD.bar(labels_d, de_summary[n_col], color="#2166ac", edgecolor="black")
    axD.set_ylabel("FDR-significant genes", fontsize=6)
    axD.tick_params(axis="x", labelsize=4.2, rotation=90)
    axD.tick_params(axis="y", labelsize=5)
    axD.set_title("D", fontsize=9, fontweight="bold", loc="left")

# ============================================================
# PANEL E
# ============================================================
example_files = [
    ("outputs_TG_final/DEGs_Female_Mutant_BIMminusPF_annotated.csv", "Female Pitx2egl1/egl1, BIM-PF GEL", axE1),
    ("outputs_TG_final/DEGs_Male_Healthy_BIMplusBAK_annotated.csv", "Male Pitx2+/+, BIM+BAK", axE2),
]
for path, label, ax_e in example_files:
    res = safe_load(path, f"E - {label}")
    if res is not None:
        lfc_col = "log2FoldChange" if "log2FoldChange" in res.columns else res.columns[1]
        padj_col = "padj" if "padj" in res.columns else [c for c in res.columns if "padj" in c.lower()][0]
        res["neglog10padj"] = -np.log10(res[padj_col].clip(lower=1e-300))
        not_sig = res[padj_col] >= PADJ_THRESH
        sig_up = (res[padj_col] < PADJ_THRESH) & (res[lfc_col] > 0)
        sig_down = (res[padj_col] < PADJ_THRESH) & (res[lfc_col] < 0)
        ax_e.scatter(res.loc[not_sig, lfc_col], res.loc[not_sig, "neglog10padj"],
                     s=2.5, color="lightgray", alpha=0.5, linewidths=0)
        ax_e.scatter(res.loc[sig_up, lfc_col], res.loc[sig_up, "neglog10padj"],
                     s=4, color="#c0392b", alpha=0.8, linewidths=0, label="Up")
        ax_e.scatter(res.loc[sig_down, lfc_col], res.loc[sig_down, "neglog10padj"],
                     s=4, color="#2166ac", alpha=0.8, linewidths=0, label="Down")
        ax_e.axhline(-np.log10(PADJ_THRESH), color="black", linewidth=0.5, linestyle="--")
        ax_e.axvline(LFC_THRESH, color="black", linewidth=0.5, linestyle="--")
        ax_e.axvline(-LFC_THRESH, color="black", linewidth=0.5, linestyle="--")
        n_up, n_down = sig_up.sum(), sig_down.sum()
        ax_e.set_title(f"{label}\n({n_up} up, {n_down} down)", fontsize=6)
        ax_e.set_xlabel("log2FC", fontsize=5.5)
        ax_e.set_ylabel("-log10 padj", fontsize=5.5)
        ax_e.tick_params(axis="both", labelsize=5)
        ax_e.legend(fontsize=4.5, loc="upper right")
axE1.text(-0.18, 1.08, "E", transform=axE1.transAxes, fontsize=9, fontweight="bold")

# ============================================================
# PANEL F
# ============================================================
interaction_summary = safe_load("outputs_TG_final/interaction_summary.csv", "F - interaction summary")
if interaction_summary is not None:
    term_col = interaction_summary.columns[0]
    n_col_f = [c for c in interaction_summary.columns if "sig" in c.lower()]
    n_col_f = n_col_f[0] if n_col_f else interaction_summary.select_dtypes(include="number").columns[0]
    axF.bar(interaction_summary[term_col], interaction_summary[n_col_f], color="#3b6fb6", edgecolor="black")
    axF.set_ylabel("Significant genes", fontsize=6)
    axF.tick_params(axis="x", labelsize=4.2, rotation=90)
    axF.tick_params(axis="y", labelsize=5)
    axF.set_title("F", fontsize=9, fontweight="bold", loc="left")

# ============================================================
# PANEL G
# ============================================================
if os.path.exists("outputs_TG_final/Supplementary_Table_TG_treatment_KEGG.xlsx"):
    xl = pd.ExcelFile("outputs_TG_final/Supplementary_Table_TG_treatment_KEGG.xlsx")
    ora_sheets = [s for s in xl.sheet_names if s.startswith("ORA_") and "Interaction" not in s]
    all_pathways = []
    for sheet in ora_sheets:
        df_sheet = pd.read_excel(xl, sheet_name=sheet)
        padj_col_g = [c for c in df_sheet.columns if "adjust" in c.lower() or "fdr" in c.lower()]
        if not padj_col_g:
            continue
        sig = df_sheet[df_sheet[padj_col_g[0]] < 0.05].copy()
        sig["source_sheet"] = sheet
        sig["neglog10FDR"] = -np.log10(sig[padj_col_g[0]].clip(lower=1e-300))
        term_col_g = "Term" if "Term" in sig.columns else sig.columns[0]
        sig["Term_display"] = sig[term_col_g]
        all_pathways.append(sig[["Term_display", "neglog10FDR", "source_sheet"]])
    if all_pathways:
        combined_g = pd.concat(all_pathways).sort_values("neglog10FDR", ascending=False).head(12)
        combined_g = combined_g.iloc[::-1]
        colors_g = ["#c0392b" if "BIM-PF" in s else "#2166ac" for s in combined_g["source_sheet"]]
        axG.barh(combined_g["Term_display"], combined_g["neglog10FDR"], color=colors_g, edgecolor="black")
        axG.set_xlabel("-log10 FDR", fontsize=6)
        axG.tick_params(axis="y", labelsize=4.5)
        axG.tick_params(axis="x", labelsize=5)
axG.set_title("G", fontsize=9, fontweight="bold", loc="left")

plt.tight_layout()
plt.savefig("outputs_TG_final/figures/FigureS4_TG_exploratory_A4.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_TG_final/figures/FigureS4_TG_exploratory_A4.pdf")