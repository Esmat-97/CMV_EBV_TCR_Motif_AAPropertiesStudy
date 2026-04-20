import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1️⃣ قراءة الملف (CSV أو Excel)
df = pd.read_csv("VDJdb_CMV_TRA_TRB_Input.csv" , sep="\t")

print(df.head())
print(df.columns) # غيّر الاسم حسب الملف اللي عندك








# استخراج الموتيفات
motif_F = df["CDR3b"].str.extract(r'^(CASS\w+)(F+)$')
motif_Y = df["CDR3b"].str.extract(r'^(CASS\w+)(Y)$')
motif_E = df["CDR3b"].str.extract(r'^(CASS\w+)(E)$')

# حساب التكرار
count_F = motif_F.dropna().shape[0]
count_Y = motif_Y.dropna().shape[0]
count_E = motif_E.dropna().shape[0]
total = df["CDR3b"].shape[0]

# تجميع النتائج في DataFrame
motif_summary = pd.DataFrame({
    "Motif": ["CASS...F", "CASS...Y", "CASS...E"],
    "Count": [count_F, count_Y, count_E],
    "Percentage": [count_F/total*100, count_Y/total*100, count_E/total*100]
})

print(motif_summary)

# رسم barplot للنسب
plt.figure(figsize=(6,4))
sns.barplot(x="Motif", y="Percentage", data=motif_summary, palette="Set2")
plt.title("Motif Distribution (CASS endings)")
plt.ylabel("Percentage (%)")
plt.savefig("CMV image/Motif Distribution (CASS endings)in CMV-specific TCRs.png")
plt.show()













# تصنيف الموتيفات
def classify_motif(seq):
    if seq.startswith("CASS"):
        if seq.endswith("F"):
            return "CASS...F"
        elif seq.endswith("Y"):
            return "CASS...Y"
        elif seq.endswith("E"):
            return "CASS...E"
    return "Other"

df["Motif"] = df["CDR3b"].apply(classify_motif)
pivot_table = df.pivot_table(index="TRBV", columns="Motif", values="CDR3b", aggfunc="count").fillna(0)

plt.figure(figsize=(12,6))
sns.heatmap(pivot_table, cmap="YlGnBu", linewidths=0.5)
plt.title("Heatmap of TRBV vs CDR3b Motifs")
plt.xlabel("CDR3b Motifs")
plt.ylabel("TRBV Gene")
plt.savefig("CMV image/Heatmap of TRBV vs CDR3b Motifs in CMV-specific TCRs.png")
plt.show()






df["Motif"] = df["CDR3b"].apply(classify_motif)

# Boxplot للأطوال حسب TRBV
plt.figure(figsize=(12,6))
sns.boxplot(x="Motif", y="TRBV", data=df, palette="Set3")
plt.xticks(rotation=45)
plt.title("Motif Distribution by TRBV Gene")
plt.ylabel("Motif")
plt.xlabel("TRBV Gene")
plt.savefig("CMV image/Motif Distribution by TRBV Gene (box Plot) in CMV-specific TCRs.png")
plt.show()




df["Motif"] = df["CDR3b"].apply(classify_motif)

# Violin plot للأطوال حسب TRBV
plt.figure(figsize=(12,6))
sns.violinplot(x="Motif", y="TRBV", data=df, palette="Set3", inner="box")
plt.xticks(rotation=45)
plt.title("Motif Distribution by TRBV Gene (Violin Plot)")
plt.ylabel("Motif")
plt.xlabel("TRBV Gene")
plt.savefig("CMV image/Motif Distribution by TRBV Gene (Violin Plot) in CMV-specific TCRs.png")
plt.show()





from scipy.stats import chi2_contingency

contingency = pd.crosstab(df["TRBV"], df["Motif"])
chi2, p, dof, ex = chi2_contingency(contingency)

print("p-value:", p)

import numpy as np

n = contingency.sum().sum()
cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape)-1)))

print("Cramér's V:", cramers_v)






