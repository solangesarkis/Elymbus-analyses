import os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt

os.makedirs("outputs_full/figures", exist_ok=True)

# ============================================================
# DONNEES : les 7 coefficients d'interaction du modele factoriel complet
# ============================================================
interaction_terms = [
    ("sex:genotype", 21, False),
    ("sex:treatment\n[BIM+BAK]", 27, False),
    ("sex:treatment\n[BIM-PF GEL]", 3, False),
    ("genotype:treatment\n[BIM+BAK]", 1, False),
    ("genotype:treatment\n[BIM-PF GEL]", 0, False),
    ("sex:genotype:treatment\n[BIM+BAK]", 2, True),
    ("sex:genotype:treatment\n[BIM-PF GEL]", 123, True),
]

labels = [t[0] for t in interaction_terms]
values = [t[1] for t in interaction_terms]
highlight = [t[2] for t in interaction_terms]

colors = ["#c9184a" if h else "#8da0cb" for h in highlight]
edge_colors = ["black" if h else "#555555" for h in highlight]
edge_widths = [1.8 if h else 0.6 for h in highlight]

# ============================================================
# FIGURE
# ============================================================
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(labels, values, color=colors, edgecolor=edge_colors, linewidth=edge_widths)

for bar, val, h in zip(bars, values, highlight):
    fontweight = "bold" if h else "normal"
    fontsize = 12 if h else 9
    ax.text(bar.get_x() + bar.get_width() / 2, val + 1.5, str(val),
            ha="center", va="bottom", fontsize=fontsize, fontweight=fontweight)

ax.set_ylabel("Number of FDR-significant genes (padj < 0.05)", fontsize=11)
ax.set_title("Two-way and three-way interaction coefficients\nof the full sex \u00d7 genotype \u00d7 treatment factorial model",
              fontsize=12)
ax.tick_params(axis="x", labelsize=8, rotation=20)

# Legende expliquant le highlight
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#c9184a", edgecolor="black", linewidth=1.5,
          label="Three-way (sex\u00d7genotype\u00d7treatment) coefficients"),
    Patch(facecolor="#8da0cb", edgecolor="#555555",
          label="Other two-way / three-way coefficients"),
]
ax.legend(handles=legend_elements, loc="upper left", fontsize=9)

plt.tight_layout()
plt.savefig("outputs_full/figures/Figure_interaction_coefficients.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_full/figures/Figure_interaction_coefficients.pdf")
