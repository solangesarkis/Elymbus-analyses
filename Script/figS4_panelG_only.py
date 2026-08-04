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

os.makedirs("outputs_TG_final/figures", exist_ok=True)

def read_csv_auto(path, **kwargs):
    try:
        return pd.read_csv(path, sep=None, engine="python", encoding="utf-8", **kwargs)
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=None, engine="python", encoding="latin-1", **kwargs)

fig, ax = plt.subplots(figsize=(4.5, 3.2))

if os.path.exists("outputs_TG_final/Supplementary_Table_TG_treatment_KEGG.xlsx"):
    xl = pd.ExcelFile("outputs_TG_final/Supplementary_Table_TG_treatment_KEGG.xlsx")
    ora_sheets = [s for s in xl.sheet_names if s.startswith("ORA_") and "Interaction" not in s]
    print(f"Feuilles ORA utilisees: {ora_sheets}")

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

        def get_color(sheet):
            if "BIM-PF" in sheet:
                return "#f4a582" if "EXPLOR" in sheet else "#d62728"
            return "#2166ac"

        colors_g = [get_color(s) for s in combined_g["source_sheet"]]
        ax.barh(combined_g["Term_display"], combined_g["neglog10FDR"],
                color=colors_g, edgecolor="black", linewidth=0.4)
        ax.set_xlabel("-log10 FDR", fontsize=6)
        ax.tick_params(axis="y", labelsize=5)
        ax.tick_params(axis="x", labelsize=5)
    else:
        ax.text(0.5, 0.5, "No significant pathways", ha="center", va="center", transform=ax.transAxes, fontsize=6)
else:
    ax.text(0.5, 0.5, "File not found", ha="center", va="center", transform=ax.transAxes, fontsize=6)

ax.set_title("G", fontsize=8, fontweight="bold", loc="left")

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#2166ac", edgecolor="black", label="BIM+BAK"),
    Patch(facecolor="#d62728", edgecolor="black", label="BIM-PF GEL"),
    Patch(facecolor="#f4a582", edgecolor="black", label="Exploratory (n=2)"),
]
ax.legend(handles=legend_elements, fontsize=4.2, loc="lower right", frameon=False)

plt.tight_layout()
plt.savefig("outputs_TG_final/figures/FigureS4_panelG_compact_A4.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_TG_final/figures/FigureS4_panelG_compact_A4.pdf")