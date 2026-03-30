# Yam GWAS (`gwas_yam.R`)

GAPIT-based association on a **numeric genotype matrix** plus map and phenotypes, with in-script QC (MAF, missingness, HWE) and post-GWAS summaries (Bonferroni, regions, tables).

**Author:** Babatunde Afeez Olabuntu

## Run

```bash
cd GWAS
Rscript gwas_yam.R
```

## Packages

| Package | Use |
|---------|-----|
| **remotes** | Install GAPIT from GitHub |
| **GAPIT** | GWAS, kinship, PCA, Manhattan/QQ |
| **utils** | Reading/writing tables and CSV |
| **stats** | SD, HWE chi-square p-values in QC |

After `library(remotes)` you can run `install_github("jiabowang/GAPIT")` once if the script fails to install it.

## What this script does vs a full server pipeline

**In this repository**, `gwas_yam.R` **does not** read VCF/BCF or call **BCFtools**, **VCFtools**, or **PLINK**. It **simulates** genotypes and phenotypes in R and writes `data/derived/geno.txt`, `map.txt`, and `pheno.txt` so the workflow runs without sharing real panels.

**On the full study**, the usual order is: raw sequence or chip data on a server → variant files (**VCF** / **BCF**) → filtering and subsetting with **BCFtools** (and often **VCFtools** for related tasks) → optional **PLINK** (or PLINK 2) for **QC** (`--maf`, `--geno`, `--mind`, `--hwe`) and export to a numeric matrix or compatible format → **then** R with **GAPIT** using the same kind of inputs this script produces after its own QC step. This README names those tools because they are the standard upstream steps; they are **not** executed by `gwas_yam.R` here.

## QC before GWAS (inside R)

The script applies **line** missingness, **SNP** MAF and missing rate, and **HWE** p-values—conceptually aligned with PLINK-style filters, but applied to the simulated matrix:

- **Lines**: drop lines with too much missing (`QC_IND_MISS_MAX`)
- **SNPs**: MAF, missing rate, HWE p-value

Filtered genotypes are written to `data/derived/geno.txt` for GAPIT.

Read **`results/qc_report.txt`** and **`results/tables/qc_per_snp_after_filters.csv`**.

## Outputs

| Path | Content |
|------|---------|
| `data/derived/pheno.txt` | Line IDs and traits (simulated here) |
| `data/derived/map.txt` | SNP, chromosome, position |
| `data/derived/geno.txt` | Taxa + dosage matrix after QC |
| `data/derived/simulation_meta.txt` | Simulation metadata (this repo only) |
| `results/qc_report.txt` | QC summary |
| `results/gapit/` | GAPIT Manhattan/QQ, GWAS CSV, PCA files |
| `results/tables/qc_per_snp_after_filters.csv` | Per-SNP QC after filters |
| `results/tables/significant_snps_all_models.csv` | Combined significant hits |
| `results/tables/regions_*.csv` | LD-style regions per result file |
| `results/STUDY_NOTES.txt` | Run notes and per-model significance lines |

