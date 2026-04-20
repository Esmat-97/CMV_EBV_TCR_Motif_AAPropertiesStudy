import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")
df = df.dropna(subset=["CDR3b","TRBV"])

def classify_property(seq):
    if not isinstance(seq, str):
        return "Other"
    aa = seq[-3]

    if aa in ["F","Y","W"]:
        return "Aromatic"
    elif aa in ["D","E"]:
        return "Negative"
    elif aa in ["K","R","H"]:
        return "Positive"
    elif aa in ["A","V","L","I","M"]:
        return "Hydrophobic"
    else:
        return "Other"

df["Property"] = df["CDR3b"].apply(classify_property)

# مهم جدًا
print(df["Property"].value_counts())

plt.figure(figsize=(6,4))
sns.countplot(data=df, x="Property")
plt.title("Physicochemical Properties of CDR3 Endings")
plt.show()

# طول CDR3
df["CDR3_length"] = df["CDR3b"].str.len()

plt.figure(figsize=(6,4))
sns.boxplot(x="Property", y="CDR3_length", data=df)
plt.title("CDR3 Length by Property")
plt.show()










# تحميل baseline
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")

def classify_property(seq):
    if not isinstance(seq, str):
        return "Other"
    aa = seq[-3]

    if aa in ["F","Y","W"]:
        return "Aromatic"
    elif aa in ["D","E"]:
        return "Negative"
    elif aa in ["K","R","H"]:
        return "Positive"
    elif aa in ["A","V","L","I","M"]:
        return "Hydrophobic"
    else:
        return "Other"

df["Property"] = df["CDR3b"].apply(classify_property)
df_human["Property"] = df_human["CDR3b"].apply(classify_property)

EBV_freq = df["Property"].value_counts(normalize=True)
human_freq = df_human["Property"].value_counts(normalize=True)

enrichment = pd.DataFrame({
    "EBV": EBV_freq,
    "Human": human_freq
}).fillna(0)

import numpy as np
enrichment["log2_enrichment"] = np.log2((enrichment["EBV"]+1e-9)/(enrichment["Human"]+1e-9))



print(enrichment)



enrichment = enrichment.sort_values("log2_enrichment", ascending=False)

print("\n=== Motif Enrichment (EBV vs Human) ===")
print(enrichment)

# رسم barplot
plt.figure(figsize=(8,5))
sns.barplot(x=enrichment.index, y=enrichment["log2_enrichment"], palette="Set2")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("Motif")
plt.title("Physicochemical Enrichment of CDR3 Endings in EBV vs Human Baseline")
plt.tight_layout()
plt.show()


