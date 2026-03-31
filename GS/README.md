# Trait — genomic selection (R) and ML/DL (Python)

This folder contains **genomic prediction** workflows for a trait, using the same genotype and phenotype inputs in **R** (mixed-model and Bayesian GS) and **Python** (classical ML and neural networks).

## Data inputs

Scripts expect these CSV files in the project root (or adjust `setwd()` / paths):

| File | Role |
|------|------|
| `dosage_boild_converted.csv` | SNP dosages (0/1/2); genotype ID in a `geno` column or as row names |
| `index_boild_mapping.csv` | Phenotypes: `geno` (line ID) and `MGIDI` (BLUP / target trait) |

R scripts often scale `MGIDI`; the notebook renames `MGIDI` → `BLUP` for consistency.

## R models

All R scripts use **10-fold cross-validation**: test folds are masked, models are fit on the rest, and predictions are evaluated with correlation **r**, **R²** (coefficient of determination), **RMSE**, and **MAE**.

### Dependencies

- **rrBLUP** — additive relationship matrix `A.mat()`, GBLUP-style fitting (`mixed.solve` / kinship BLUP).
- **tidyverse** — data handling (many scripts).
- **sommer** — mixed models with multiple variance components (some GBLUP and RKHS setups); see `GBLUP_Comprehensive_Models.R`.
- **BGLR** — Bayesian alphabet models (BayesA/B/C/Cπ/R).
- **randomForest** — where used in the comprehensive script.

Install from CRAN as needed, e.g. `install.packages(c("rrBLUP", "tidyverse", "sommer", "BGLR", "randomForest"))`.

### Script map

| Script | Method | Typical output under `Results/` |
|--------|--------|----------------------------------|
| `additive.R` | Additive GBLUP | `additivegblup_10folds.csv`, variance components CSV |
| `additive_dominance.R` | Additive + dominance GBLUP | `10foldcv_additive_dominance_gblup.csv` |
| `Additive_epistatic.R` | Additive + epistatic GBLUP | `10foldcv_additive_epistatic_gblup.csv` |
| `Dominance.R` | Dominance GBLUP | `10foldcv_dominance_gblup.csv` |
| `Additive_GXE.R` | Additive GBLUP with G×E | `10foldcv_additive_GXE_gblup.csv` |
| `RKHS.R` | Reproducing kernel Hilbert space | `10foldcv_rkhs.csv` |
| `BayesA.R`, `BayesB.R`, `BayesC.R`, `BayesCpi.R`, `BayesR.R` | Bayesian GS (BGLR) | `10foldcv_bayesa.csv`, … `bayesr.csv` |
| `GBLUP_Comprehensive_Models.R` | Several GBLUP / sommer / RF variants in one place | (per script outputs; see file headers) |

**Note:** `compare_all_methods.R` expects additive GBLUP metrics in `Results/10foldcv_additivegblup.csv`. The standalone `additive.R` currently writes `additivegblup_10folds.csv`. Rename or copy to match before running the comparison, or adjust `compare_all_methods.R`.

### Aggregating all methods

After R and Python runs have produced the per-method CSVs in `Results/`, run:

```r
source("compare_all_methods.R")
```

This reads GBLUP, Bayesian, RKHS, ML, and DL result files, computes mean **r**, **R²**, **RMSE**, **MAE** per method, sorts by **R²**, and writes `Results/comparison_all_methods.csv`.

Update `setwd()` at the top of each script if your project lives elsewhere.

---

## `ML_and_DL.ipynb` (Python)

Jupyter notebook for **machine learning** and **deep learning** on the same `dosage_boild_converted.csv` and `index_boild_mapping.csv` data.

### Environment

Typical stack:

- **pandas**, **numpy**, **matplotlib**, **seaborn**
- **scikit-learn** — preprocessing, `KFold`, metrics (`r2_score`, MSE, MAE)
- **xgboost**; optionally **lightgbm**, **lazypredict** (`LazyRegressor`)
- **PyTorch** — feedforward NN and CNN
- **RDKit** — imported in early cells (used if cheminformatics features are explored; core GS cells use genotype matrices)

Install matching versions for your Python (e.g. conda or pip). GPU is optional for the small networks shown.

### What the notebook does (flow)

1. **Load and align** genotype and phenotype tables; check ID alignment and duplicates.
2. **Preprocess** — scaling options, merges, optional **feature selection** (e.g. variance filters, univariate tests) and MAF-style filtering in later diagnostic sections.
3. **Baseline / exploration** — includes **LazyRegressor**-style screening of many sklearn regressors.
4. **Classical ML with 10-fold CV** — e.g. SVM, Random Forest, XGBoost, LightGBM, Gradient Boosting, Extra Trees, k-NN, AdaBoost; metrics aligned with the R pipeline (**r**, **R²**, RMSE, MAE). Results are written under `Results/` as `svm_results.csv`, `rf_results.csv`, `xgb_results.csv`, etc.
5. **Deep learning** — **feedforward (FWNN)** and **CNN** sections on marker matrices; training loops with PyTorch.
6. **Diagnostics** — sections on negative **R²**, phenotype scaling, heritability-style checks, Ridge/Lasso/ElasticNet, and a “genomic prediction style” pipeline (e.g. MAF filtering, ridge-style regression) to interpret difficult prediction scenarios.




