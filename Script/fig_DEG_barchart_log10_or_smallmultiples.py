import os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("outputs_full/figures", exist_ok=True)

# ============================================================
# DONNEES (identiques a la version precedente)
# ============================================================
comparisons = [
    ("Female Pitx2+/+\nBIM+BAK vs Untreated", 516, -872),
    ("Female Pitx2+/+\nBIM-PF GEL vs Untreated", 2, -13),
    ("Male Pitx2+/+\nBIM+BAK vs Untreated", 112, -95),
    ("Male Pitx2+/+\nBIM-PF GEL vs Untreated", 80, -62),
    ("Female Pitx2egl1/egl1\nBIM+BAK vs Untreated", 58, -55),
    ("Female Pitx2egl1/egl1\nBIM-PF GEL vs Untreated", 11, -1),
    ("Male Pitx2egl1/egl1\nBIM+BAK vs Untreated", 3, -11),
    ("Male Pitx2egl1/egl1\nBIM-PF GEL vs Untreated", 0, -1),
]

labels = [c[0] for c in comparisons]
up_vals = [c[1] for c in comparisons]
down_vals = [c[2] for c in comparisons]
x = range(len(labels))

# ============================================================
# OPTION A : UN SEUL PANNEAU, ECHELLE log10(count+1), diverging bar chart
# ============================================================
def signed_log10(v):
    """log10(|v|+1), en gardant le signe d'origine."""
    return np.sign(v) * np.log10(np.abs(v) + 1)

up_log = [signed_log10(v) for v in up_vals]
down_log = [signed_log10(v) for v in down_vals]

fig, ax = plt.subplots(figsize=(11, 6))
ax.bar(x, up_log, color="#c0392b", edgecolor="black", linewidth=0.5, label="Upregulated")
ax.bar(x, down_log, color="#2166ac", edgecolor="black", linewidth=0.5, label="Downregulated")
ax.axhline(0, color="black", linewidth=0.8)

# Annotation des vraies valeurs (pas des valeurs log) au-dessus/en-dessous de chaque barre
for i, (u, d) in enumerate(zip(up_vals, down_vals)):
    u_log = signed_log10(u)
    d_log = signed_log10(d)
    ax.text(i, u_log + 0.05, str(u), ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.text(i, d_log - 0.05, str(abs(d)), ha="center", va="top", fontsize=8, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
ax.set_ylabel("log$_{10}$(DEG count + 1), signed by direction", fontsize=10)
ax.set_title("Differentially expressed genes per comparison (treatment vs untreated)\n"
              "log$_{10}$(count+1) scale, no axis break", fontsize=11)
ax.legend(loc="upper right", fontsize=9)

plt.tight_layout()
plt.savefig("outputs_full/figures/Figure_DEG_barchart_log10.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_full/figures/Figure_DEG_barchart_log10.pdf")

# ============================================================
# OPTION B : PETITS MULTIPLES ALIGNES (une sous-figure par comparaison,
# meme echelle Y fixe et identique partout -- alternative au log10)
# ============================================================
fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey=False)
axes = axes.flatten()

y_max = max(max(up_vals), max(abs(v) for v in down_vals)) * 1.15

for ax, label, u, d in zip(axes, labels, up_vals, down_vals):
    ax.bar(["Up", "Down"], [u, abs(d)], color=["#c0392b", "#2166ac"], edgecolor="black")
    ax.text(0, u, str(u), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.text(1, abs(d), str(abs(d)), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_ylim(0, y_max)
    ax.set_title(label, fontsize=8)
    ax.tick_params(axis="both", labelsize=8)

fig.suptitle("Differentially expressed genes per comparison (treatment vs untreated)\n"
              "Aligned small multiples, common linear y-axis", fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig("outputs_full/figures/Figure_DEG_barchart_small_multiples.pdf", bbox_inches="tight")
plt.close()
print("Sauvegarde: outputs_full/figures/Figure_DEG_barchart_small_multiples.pdf")
