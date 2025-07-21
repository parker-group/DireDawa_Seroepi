# High Seroprevalence of Antibodies to Dengue, Chikungunya, and Zika Viruses in Dire Dawa, Ethiopia: A Cross-Sectional Survey in 2024

This repository contains the full analytical workflow and data used in our cross-sectional serological survey conducted in Dire Dawa, Ethiopia. The study examines IgG and IgM responses to dengue virus (DENV), chikungunya virus (CHIKV), and Zika virus (ZIKV) in 2024.

---

## Repository Structure

```
DireDawa_Seroepi/
├── data/
│   ├── DireDawa_SeroEpi24.csv         # Primary dataset (Dire Dawa)
│   └── AddisAbaba_SeroEpi24.csv       # Comparison dataset (Addis Ababa)
├── docs/
│   └── data_dictionary.md             # Description of variables used in the dataset
├── scripts/
│   ├── OverallSeroprevEstimates.r     # Overall IgG and IgM seroprevalence estimates with bootstrapped CIs
│   ├── AgeSpecificSeroDireDawa.r      # Age-specific seropositivity and plots (Dire Dawa)
│   ├── AgeSpecificSero_Addis.r        # Same analysis for Addis Ababa data
│   ├── GEE_AgeSpecificSero.r          # GEE with splines to model seropositivity by age
│   └── CrossReactivityICC.r           # Mixed models, correlation matrix, and ICC clustering estimates
└── README.md
```

---

## Scripts and Workflow

### 1. [`scripts/OverallSeroprevEstimates.r`](scripts/OverallSeroprevEstimates.r)
Estimates overall IgG and IgM seroprevalence for DENV, CHIKV, and ZIKV. Includes bootstrap confidence intervals based on both individual and household-level resampling.

### 2. [`scripts/AgeSpecificSeroDireDawa.r`](scripts/AgeSpecificSeroDireDawa.r)
Calculates seroprevalence by age group in Dire Dawa and generates bar plots for IgG and IgM.

### 3. [`scripts/AgeSpecificSero_Addis.r`](scripts/AgeSpecificSero_Addis.r)
Performs the same age-specific analysis on the Addis Ababa comparison dataset.

### 4. [`scripts/GEE_AgeSpecificSero.r`](scripts/GEE_AgeSpecificSero.r)
Fits Generalized Estimating Equation (GEE) models with natural splines to estimate age-specific seropositivity trends.

### 5. [`scripts/CrossReactivityICC.r`](scripts/CrossReactivityICC.r)
- Uses linear mixed models to evaluate cross-reactivity between DENV, ZIKV, and CHIKV IgG and IgM values.
- Generates a Spearman correlation heatmap across all antibody pairs.
- Calculates intra-class correlation coefficients (ICCs) to evaluate household clustering of antibody responses.

---

## Data

- **Primary dataset:**  
  [`data/DireDawa_SeroEpi24.csv`](data/DireDawa_SeroEpi24.csv) — Individual-level serological and demographic data for Dire Dawa, Ethiopia.

- **Comparison dataset:**  
  [`data/AddisAbaba_Seroepi24.csv`](data/AddisAbaba_Seroepi24.csv) — Similar dataset from Addis Ababa for comparative analysis.

- **Data dictionary:**  
  Variable definitions and metadata are provided in [`docs/data_dictionary.md`](docs/data_dictionary.md).

---

## Citation

If you use this repository, please cite:

> **Parker DM, et al.** High Seroprevalence of Antibodies to Dengue, Chikungunya, and Zika Viruses in Dire Dawa, Ethiopia: A Cross-Sectional Survey in 2024. Preprint DOI: 10.1101/2024.09.19.24313917


