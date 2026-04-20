import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



# -------------------------------
# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_EBV = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# دالة لتصنيف الموتيفات
def classify_motif(seq):
    if seq.startswith("CASS"):
        if seq.endswith("F"):
            return "CASS...F"
        elif seq.endswith("Y"):
            return "CASS...Y"
        elif seq.endswith("E"):
            return "CASS...E"
    return "Other"

df_human["Motif"] = df_human["CDR3b"].apply(classify_motif)
df_EBV["Motif"] = df_EBV["CDR3b"].apply(classify_motif)

# حساب frequency لكل motif
human_freq = df_human["Motif"].value_counts(normalize=True)
EBV_freq = df_EBV["Motif"].value_counts(normalize=True)

# دمج الاتنين
df_motif = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": EBV_freq}).fillna(0)

# حساب log2 enrichment
df_motif["log2_enrichment"] = np.log2((df_motif["EBV_freq"]+1e-9) / (df_motif["Human_freq"]+1e-9))

# ترتيب
df_motif = df_motif.sort_values("log2_enrichment", ascending=False)

print("\n=== Motif Enrichment (EBV vs Human) ===")
print(df_motif)

# رسم barplot
plt.figure(figsize=(8,5))
sns.barplot(x=df_motif.index, y=df_motif["log2_enrichment"], palette="Set2")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("Motif")
plt.title("EBV Motif Enrichment in EBV relative to Human baseline")
plt.tight_layout()
plt.show()











