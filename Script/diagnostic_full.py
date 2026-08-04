import pandas as pd

files = {
    ("Female", "Healthy", "Untreated"): "data/Cornea Read counts Female Healthy Untreated 1.xlsx",
    ("Female", "Mutant",  "Untreated"): "data/Cornea Read counts Female Mutant Untreated 1.xlsx",
    ("Male",   "Healthy", "Untreated"): "data/Cornea Read counts Male Healthy Untreated 1.xlsx",
    ("Male",   "Mutant",  "Untreated"): "data/Cornea Read counts Male Mutant Untreated 1.xlsx",
    ("Female", "Healthy", "BIM+BAK"):   "data/Cornea Read counts Female Control BIM+BAK.xlsx",
    ("Female", "Mutant",  "BIM+BAK"):   "data/Cornea Read counts Female Mutant BIM+BAK.xlsx",
    ("Male",   "Healthy", "BIM+BAK"):   "data/Cornea Read counts Male Control BIM+BAK.xlsx",
    ("Male",   "Mutant",  "BIM+BAK"):   "data/Cornea Read counts Male Mutant BIM+BAK.xlsx",
    ("Female", "Healthy", "BIM-PF"):    "data/Cornea Read counts Female Control BIM-PF.xlsx",
    ("Female", "Mutant",  "BIM-PF"):    "data/Cornea Read counts Female Mutant BIM-PF.xlsx",
    ("Male",   "Healthy", "BIM-PF"):    "data/Cornea Read counts Male Control BIM-PF.xlsx",
    ("Male",   "Mutant",  "BIM-PF"):    "data/Cornea Read counts Male Mutant BIM-PF.xlsx",
}

gene_sets = {}
n_cols = {}

for key, path in files.items():
    df = pd.read_excel(path)
    label = "_".join(key)
    print(f"=== {label} ===")
    print(f"Lignes (genes): {df.shape[0]}")
    id_cols = [c for c in ("Gene ID", "Gene Symbol", "Type") if c in df.columns]
    count_cols = [c for c in df.columns if c not in id_cols]
    print(f"Colonnes de comptage ({len(count_cols)}): {count_cols}")
    gene_sets[label] = set(df["Gene ID"])
    n_cols[label] = len(count_cols)

print("\n=== Verification: tous les fichiers ont-ils le meme nombre de genes ? ===")
sizes = {k: len(v) for k, v in gene_sets.items()}
print(sizes)

print("\n=== Verification: chevauchement des genes (vs premier fichier) ===")
ref_label = list(gene_sets.keys())[0]
ref_set = gene_sets[ref_label]
for label, s in gene_sets.items():
    common = len(ref_set & s)
    print(f"{ref_label} vs {label}: {common} genes communs / {len(ref_set)} et {len(s)}")

print("\n=== Nombre total d'echantillons attendu ===")
print(f"Somme des colonnes de comptage: {sum(n_cols.values())} (attendu: 60 si 5 replicats x 12 groupes)")
