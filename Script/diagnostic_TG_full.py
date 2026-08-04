import pandas as pd

files = {
    ("Female", "Healthy", "Untreated"): "data/Female Control NT TG.xlsx",
    ("Female", "Mutant",  "Untreated"): "data/Female Mutant NT TG.xlsx",
    ("Male",   "Healthy", "Untreated"): "data/Male Control NT TG.xlsx",
    ("Male",   "Mutant",  "Untreated"): "data/Male Mutant NT TG.xlsx",
    ("Female", "Healthy", "BIM+BAK"):   "data/Female Control BIM+BAK TG.xlsx",
    ("Female", "Mutant",  "BIM+BAK"):   "data/Female Mutant BIM+BAK TG.xlsx",
    ("Male",   "Healthy", "BIM+BAK"):   "data/Male Control BIM+BAK TG.xlsx",
    ("Male",   "Mutant",  "BIM+BAK"):   "data/Male Mutant BIM+BAK TG.xlsx",
    ("Female", "Healthy", "BIM-PF"):    "data/Female Control BIM-PF TG.xlsx",
    ("Female", "Mutant",  "BIM-PF"):    "data/Female Mutant BIM-PF TG.xlsx",
    ("Male",   "Healthy", "BIM-PF"):    "data/Male Control BIM-PF TG.xlsx",
    ("Male",   "Mutant",  "BIM-PF"):    "data/Male Mutant BIM-PF TG.xlsx",
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

print("\n=== Verification: meme nombre de genes partout ? ===")
sizes = {k: len(v) for k, v in gene_sets.items()}
print(sizes)

print("\n=== Verification: chevauchement (vs premier fichier) ===")
ref_label = list(gene_sets.keys())[0]
ref_set = gene_sets[ref_label]
for label, s in gene_sets.items():
    common = len(ref_set & s)
    print(f"{ref_label} vs {label}: {common} genes communs / {len(ref_set)} et {len(s)}")

print("\n=== Nombre total d'echantillons attendu ===")
print(f"Somme des colonnes de comptage: {sum(n_cols.values())} (attendu: 60 si 5 replicats x 12 groupes)")
