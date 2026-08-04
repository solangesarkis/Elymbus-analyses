import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

os.makedirs("outputs_full", exist_ok=True)

counts_for_deseq = pd.read_csv("outputs_full/counts_filtered_for_deseq.csv", index_col=0)
print(f"Matrice: {counts_for_deseq.shape[0]} echantillons x {counts_for_deseq.shape[1]} genes")

meta_rows = []
for sample in counts_for_deseq.index:
    parts = sample.split("_")
    sex, genotype, treatment = parts[0], parts[1], parts[2]
    meta_rows.append({"sample": sample, "sex": sex, "genotype": genotype, "treatment": treatment})
metadata = pd.DataFrame(meta_rows).set_index("sample").loc[counts_for_deseq.index]
print(metadata.groupby(["sex", "genotype", "treatment"]).size())

inference = DefaultInference(n_cpus=4)

# ============================================================
# PARTIE A : COMPARAISONS PAR SOUS-GROUPE
# (pour chaque combinaison sexe x genotype : chaque traitement vs Untreated)
# ============================================================
print("\n\n########## PARTIE A: comparaisons par sous-groupe ##########")
summary_rows = []

for sex in ["Female", "Male"]:
    for genotype in ["Healthy", "Mutant"]:
        subset_mask = (metadata["sex"] == sex) & (metadata["genotype"] == genotype)
        sub_counts = counts_for_deseq.loc[subset_mask]
        sub_meta = metadata.loc[subset_mask]

        print(f"\n=== {sex} - {genotype}: modele ~treatment ===")
        dds_sub = DeseqDataSet(
            counts=sub_counts,
            metadata=sub_meta,
            design="~treatment",
            refit_cooks=True,
            inference=inference,
        )
        dds_sub.deseq2()

        for treat in ["BIM+BAK", "BIM-PF"]:
            stat = DeseqStats(dds_sub, contrast=["treatment", treat, "Untreated"], inference=inference)
            stat.summary()
            label = f"{sex}_{genotype}_{treat.replace('+','plus').replace('-','minus')}_vs_Untreated"
            stat.results_df.to_csv(f"outputs_full/DE_{label}.csv")
            n_sig = (stat.results_df["padj"] < 0.05).sum()
            n_up = ((stat.results_df["padj"] < 0.05) & (stat.results_df["log2FoldChange"] > 0)).sum()
            n_down = ((stat.results_df["padj"] < 0.05) & (stat.results_df["log2FoldChange"] < 0)).sum()
            print(f"{sex} {genotype}: {treat} vs Untreated -> {n_sig} DEGs ({n_up} up, {n_down} down)")
            summary_rows.append({"sex": sex, "genotype": genotype, "comparison": f"{treat} vs Untreated",
                                  "n_sig": n_sig, "n_up": n_up, "n_down": n_down})

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("outputs_full/DE_summary_by_subgroup.csv", index=False)
print("\n--- RESUME PARTIE A ---")
print(summary_df.to_string(index=False))

# ============================================================
# PARTIE B : MODELE FACTORIEL COMPLET (triple interaction)
# ============================================================
print("\n\n########## PARTIE B: modele factoriel complet (60 echantillons) ##########")
dds_full = DeseqDataSet(
    counts=counts_for_deseq,
    metadata=metadata,
    design="~sex + genotype + treatment + sex:genotype + sex:treatment + genotype:treatment + sex:genotype:treatment",
    refit_cooks=True,
    inference=inference,
)
dds_full.deseq2()

design_cols = list(dds_full.obsm["design_matrix"].columns)
print("\n--- Colonnes de la matrice de design (modele complet) ---")
for c in design_cols:
    print(" -", c)

print("\n--- Facteurs de taille ---")
print(dds_full.obs["size_factors"].describe())

# Extraire chaque coefficient d'interaction (2-way et 3-way) automatiquement
interaction_cols = [c for c in design_cols if ":" in c]
print(f"\n{len(interaction_cols)} colonnes d'interaction detectees.")

interaction_summary = []
for col in interaction_cols:
    idx = design_cols.index(col)
    contrast_vector = np.zeros(len(design_cols))
    contrast_vector[idx] = 1
    stat = DeseqStats(dds_full, contrast=contrast_vector, inference=inference)
    stat.summary()
    safe_name = col.replace("[", "_").replace("]", "").replace(":", "_X_").replace(".", "")
    stat.results_df.to_csv(f"outputs_full/DE_interaction_{safe_name}.csv")
    n_sig = (stat.results_df["padj"] < 0.05).sum()
    print(f"{col}: {n_sig} genes significatifs (padj<0.05)")
    interaction_summary.append({"interaction_term": col, "n_significant": n_sig})

pd.DataFrame(interaction_summary).to_csv("outputs_full/interaction_summary.csv", index=False)
print("\n=== TERMINE ===")
