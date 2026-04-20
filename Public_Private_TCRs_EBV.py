import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. قراءة البيانات
df = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# 2. تصنيف الموتيفات بالنهايات F/Y/E/Other
def classify_motif(seq):
    if isinstance(seq, str) and seq.startswith("CASS"):
        if seq.endswith("F"):
            return "CASS...F"
        elif seq.endswith("Y"):
            return "CASS...Y"
        elif seq.endswith("E"):
            return "CASS...E"
    return "Other"

df["Motif"] = df["CDR3b"].apply(classify_motif)

# 3. تقسيم Public vs Private
public = df.groupby("CDR3b").filter(lambda x: x["subject"].nunique() > 1)
private = df.groupby("CDR3b").filter(lambda x: x["subject"].nunique() == 1)

public["Category"] = "Public"
private["Category"] = "Private"

df_all = pd.concat([public, private])

# 4. جدول النسب
summary_table = df_all.groupby(["Motif","Category"])["CDR3b"].nunique().reset_index(name="Count")
summary_table["Percentage"] = summary_table.groupby("Category")["Count"].transform(lambda x: 100 * x / x.sum())

print(summary_table)

# 5. Barplot للمقارنة بين Public و Private
plt.figure(figsize=(10,6))
sns.barplot(data=summary_table, x="Motif", y="Percentage", hue="Category")
plt.ylabel("Percentage (%)")
plt.title("EBV Motif distribution in Public vs Private TCRs")
plt.legend(title="Category")
plt.show()
