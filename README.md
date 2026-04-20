# 🧬 CMV TCR Repertoire Analysis

## 📌 Overview

This project presents a computational analysis of CMV-associated T-cell receptor (TCR) repertoires using paired α/β chain data. The goal is to explore sequence-level patterns, gene associations, and physicochemical properties of CDR3 regions.

---

## 🎯 Objectives

* Analyze CDR3β motif patterns (e.g., CASS…F/Y/E)
* Compare CMV-specific repertoires with a human baseline dataset
* Identify motif enrichment and depletion patterns
* Investigate public vs private clonotypes
* Explore physicochemical properties of terminal amino acids
* Evaluate TRBV/TRAV gene associations with motif classes

---

## 🧪 Methods

### 1. Data Processing

* Input: VDJdb-derived CMV and human baseline datasets
* Cleaning: removal of invalid sequences and normalization of gene names
* Separation of TRA and TRB chains with paired analysis

### 2. Motif Analysis

* Classification of CDR3β sequences into:

  * CASS…F
  * CASS…Y
  * CASS…E
  * Other
* Frequency and percentage calculations

### 3. Enrichment Analysis

* Computation of log2 enrichment:

  log2(CMV frequency / Human frequency)

* Identification of overrepresented and underrepresented motifs

### 4. Public vs Private Clonotypes

* Public: shared across multiple subjects
* Private: unique to a single subject
* Comparative distribution analysis

### 5. Physicochemical Profiling

* Grouping terminal amino acids into:

  * Aromatic
  * Hydrophobic
  * Charged (Positive/Negative)
  * Other
* Frequency and enrichment analysis

### 6. Statistical Analysis

* Chi-square test for association between:

  * TRBV/TRAV genes and motif classes
* Effect size estimation using Cramér’s V

---

## 📊 Key Findings

* CASS…F motifs are dominant but not specifically enriched in CMV
* CMV repertoires show enrichment in hydrophobic and charged residues
* Negative residues exhibit strong relative enrichment despite low frequency
* Public TCRs show stronger motif conservation than private ones
* Gene–motif associations suggest structured repertoire organization

---

## 🧠 Interpretation

The results indicate that CMV-associated TCR repertoires are shaped not only by sequence motifs but also by underlying physicochemical constraints and gene usage patterns.

---

## ⚠️ Limitations

* Pairing of TRA/TRB chains is based on dataset annotations
* Low-frequency categories may inflate enrichment values
* Motif abstraction simplifies full sequence complexity

---

## 🛠 Tools & Libraries

* Python
* pandas
* numpy
* matplotlib
* seaborn
* scipy

---

## 📂 Output

* Bar plots (motif distribution, enrichment)
* Heatmaps (TRBV vs motif)
* Public vs private comparisons
* Physicochemical distribution plots

---

## 🚀 Future Work

* Position-specific amino acid analysis
* Entropy and diversity metrics
* Machine learning classification of TCR specificity
* Cross-pathogen comparison (CMV vs other viruses)

---

## 📬 Contact

For collaboration or discussion, feel free to connect.
