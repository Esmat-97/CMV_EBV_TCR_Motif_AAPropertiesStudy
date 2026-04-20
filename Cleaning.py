import pandas as pd
import json

# =========================
# 1. Load data
# =========================
input_file = "VDJdb_CMV_TRA_TRB.tsv"
df = pd.read_csv(input_file, sep="\t")

# =========================
# 2. Clean gene names
# =========================
def clean_gene(gene):
    if pd.isna(gene):
        return None
    return str(gene).split('*')[0].split('/')[0]

df['V'] = df['V'].apply(clean_gene)
df['J'] = df['J'].apply(clean_gene)

# =========================
# 3. Clean CDR3 (basic QC)
# =========================
df = df.dropna(subset=['CDR3'])
df = df[df['CDR3'].str.len() > 5]

# remove obvious invalid characters (optional but safe)
df = df[df['CDR3'].str.contains("^[ACDEFGHIKLMNPQRSTVWY]+$")]

# =========================
# 4. Extract subject safely
# =========================
def extract_subject(meta_str):
    try:
        meta_dict = json.loads(meta_str.replace("'", '"'))
        return meta_dict.get('subject.id') or meta_dict.get('subject.cohort') or "Unknown"
    except:
        return "Unknown"

if 'Meta' in df.columns:
    df['subject'] = df['Meta'].apply(extract_subject)
else:
    df['subject'] = "Unknown"

# =========================
# 5. Split chains cleanly (NO merging)
# =========================
beta = df[df['Gene'] == 'TRB'].copy()
alpha = df[df['Gene'] == 'TRA'].copy()

# =========================
# 6. Beta cleaning
# =========================
beta = beta[['CDR3', 'V', 'J', 'subject', 'complex.id']].copy()
beta.columns = ['CDR3b', 'TRBV', 'TRBJ', 'subject', 'complex.id']

# =========================
# 7. Alpha cleaning
# =========================
alpha = alpha[['CDR3', 'V', 'J', 'complex.id']].copy()
alpha.columns = ['CDR3a', 'TRAV', 'TRAJ', 'complex.id']

# =========================
# 8. OPTIONAL: keep both but DO NOT force pairing
# =========================
# We only attach alpha info where it naturally exists
final_df = beta.merge(alpha, on='complex.id', how='left')

# =========================
# 9. Add count column (useful for downstream stats)
# =========================
final_df['count'] = 1

# =========================
# 10. Remove exact duplicates (safe dedup only)
# =========================
final_df = final_df.drop_duplicates()

# =========================
# 11. Save
# =========================
output_file = "VDJdb_CMV_TRA_TRB_Input.csv"
final_df.to_csv(output_file, sep="\t", index=False)

print("DONE ✔ Clean dataset saved:", output_file)
print("Rows:", len(final_df))