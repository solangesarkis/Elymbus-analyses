import os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

os.makedirs("outputs_TG_final/figures", exist_ok=True)

# ============================================================
# LES 8 CONTRASTES TG
# ============================================================
comparisons = [
    ("Female\nPitx2+/+", "BIM+BAK", 0, False),
    ("Female\nPitx2+/+", "BIM-PF GEL", 0, False),
    ("Male\nPitx2+/+", "BIM+BAK", 55, False),
    ("Male\nPitx2+/+", "BIM-PF GEL", 39, True),
    ("Female\nPitx2egl1/egl1", "BIM+BAK", 0, False),
    ("Female\nPitx2egl1/egl1", "BIM-PF GEL", 16, False),
    ("Male\nPitx2egl1/egl1", "BIM+BAK", 0, False),
    ("Male\nPitx2egl1/egl1", "BIM-PF GEL", 0, False),
]

labels = [f"{c[0]}\n{c[1]}" for c in comparisons]
values = [c[2] for c in comparisons]
exploratory = [c[3] for c in comparisons]

colors = []
for treat, is_explor in zip([c[1] for c in comparisons], exploratory):
    if is_explor:
        colors.append("#f4a582")
    elif treat == "BIM+BAK":
        colors.append("#2166ac")
    else:
        colors.append("#d62728")

# ============================================================
# TAILLE COMPACTE : ~1/4 de page A4 (8.27 x 2.6 pouces environ)
# ============================================================
fig, ax = plt.subplots(figsize=(4.0, 3.0))
bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.4)

for bar, val, is_explor in zip(bars, values, exploratory):
    label_text = f"{val}*" if is_explor else str(val)
    ax.text(bar.get_x() + bar.get_width() / 2, val + 1.2, label_text,
            ha="center", va="bottom", fontsize=5.5, fontweight="bold")

ax.set_ylabel("FDR-significant genes", fontsize=6)
ax.tick_params(axis="x", labelsize=4.2, rotation=90)
ax.tick_params(axis="y", labelsize=5)
ax.set_title("D", fontsize=8, fontweight="bold", loc="left")

legend_elements = [
    Patch(facecolor="#2166ac", edgecolor="black", label="BIM+BAK"),
    Patch(facecolor="#d62728", edgecolor="black", label="BIM-PF GEL"),
    Patch(facecolor="#f4a582", edgecolor="black", label="Exploratory (n=2)"),
]
ax.legend(handles=legend_elements, fontsize=4.2, loc="upper left", frameon=False)

plt.tight_layout()
plt.savefig("outputs_TG_final/figures/FigureS4_panelD_compact_A4.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_TG_final/figures/FigureS4_panelD_compact_A4.pdf")