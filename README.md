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

1. Stability Factor
Proteins with an Instability Index > 20 (conservative threshold) are penalized.
```math
Stability_Factor = 1.0 - (max(0, Instability_Index - 20) / 100.0)

2. Adsorption Factor
Proteins with high hydrophobicity (GRAVY) are assumed to have higher loss due to surface adsorption.

코드 스니펫

Adsorption_Factor = 1.0 - (abs(GRAVY) * 0.2)

3. Effective Concentration
The final recovery is a product of the Kit Efficiency (default 0.5) and the Protein's intrinsic efficiency.

코드 스니펫

Total_Coeff = Kit_Efficiency * Stability_Factor * Adsorption_Factor
C_effective = C_theoretical_max * Total_Coeff


🚀 Installation & Usage
1. Prerequisites
You need Python 3.x and Biopython.

Bash

pip install biopython
2. How to Run
Place your sequence files (.fasta or .txt) in the same directory as the script.

Run the script:

Bash

python eprm_analyzer.py
3. Check Results
A new folder named EPRM_Results_YYYYMMDD_... will be created. Check the eprm_analysis_detail.log file inside.

📂 Input File Format
Standard FASTA format is supported.

Plaintext

>Protein_A
MKVIFLKDVKGMGKKGEIKNVADGYANNFLFKQGLAIEA
>Protein_B
MKA...
📄 Output Log Example
Plaintext

[Target: Protein_A]
  • Molecular Weight: 42.50 kDa
  • Instability Index: 45.20 (Threshold: 40, Conservative: 20)
  • GRAVY (Hydrophobicity): -0.150
  • Calculated Recovery Coeff: 0.3650 (Kit: 0.5)
  • >> Estimated Effective Conc: 1.2500 uM
  • [Exp Guide] For 20nM final conc, dilute 1:62
⚙️ Configuration
You can adjust the default experimental parameters in the __init__ method of the EPRMAnalyzer class:

Python

def __init__(self, initial_conc_um=10.0, ...):
    self.eta_kit = 0.50  # Adjust kit efficiency here
⚠️ Disclaimer
This tool provides an estimation based on sequence properties. Actual recovery rates may vary depending on temperature, buffer conditions, and handling errors. Always verify concentrations with a spectrophotometer (e.g., Nanodrop) before critical experiments.

License
This project is licensed under the MIT License.
