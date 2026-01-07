# EPRM Analyzer (Effective Protein Recovery Mass Analyzer)

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Biopython](https://img.shields.io/badge/Dependency-Biopython-green)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

**EPRM Analyzer** is a Python tool designed to predict the **Effective Protein Recovery Mass (EPRM)** after purification or labeling processes. It integrates physicochemical property analysis (via Biopython) with experimental kit efficiency to provide a realistic estimation of the final protein concentration.

It goes beyond simple dilution calculations by applying a heuristic damping algorithm based on the protein's **Instability Index** and **Hydrophobicity (GRAVY)**.

---

## 📌 Key Features

- **Batch Processing:** Automatically detects and processes all `.fasta` and `.txt` files in the directory.
- **Physicochemical Analysis:** Calculates Molecular Weight (MW), Instability Index, and GRAVY using Biopython.
- **Correction Algorithm:** Applies "Stability Factor" and "Adsorption Factor" to estimate realistic protein loss.
- **Experimental Guide:** Automatically calculates the optimal dilution factor for sensitive assays (e.g., Target 20nM for MST/SPR).
- **Auto-Logging:** Creates timestamped directories and detailed logs for every run.

---

## 🧪 Scientific Logic (Calculation Model)

The tool calculates the Effective Concentration ($C_{eff}$) using the following logic:

### 1. Stability Factor
Proteins with an Instability Index > 20 (conservative threshold) are penalized.
```math
Stability_Factor = 1.0 - (max(0, Instability_Index - 20) / 100.0)
