# Population genetics analysis (R)

R workflow for summarising SNP data, choosing cluster number *K*, and producing population-structure figures after genotype data have been prepared elsewhere.

**Author:** Babatuntude Afeez Olabuntu

---

## Role in the pipeline

This script is meant to run **after** heavy lifting on a server or cluster. Upstream steps typically use command-line genomics tools on **VCF** / **BCF** (or related) data, then export formats that R can read. Common tools in that phase include:

| Tool | Typical use in this workflow |
|------|-------------------------------|
| **[BCFtools](https://samtools.github.io/bcftools/)** | View, filter, merge, and index VCF/BCF; subset samples or regions; quality and allele-count summaries. |
| **[VCFtools](https://vcftools.github.io/)** | Filtering, missingness, allele frequency, and other VCF-based statistics (where still used in your pipeline). |
| **[SAMtools](https://www.htslib.org/)** | BAM/CRAM and general HTSlib stack; often paired with BCFtools for indexed variant files. |
| **[PLINK](https://www.cog-genomics.org/plink/)** or **PLINK 2** | VCF → PLINK binary or text; SNP filtering (e.g. MAF, missingness); **allele frequencies** (`.frq`), **Hardy–Weinberg** (`.hwe`), **long-format genotypes** (`.raw`), **FAM** (`.fam`), and **genetic distance** matrices for trees. |
| **[BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html)** | Genotype **imputation** and phasing; outputs feed back into PLINK-style files for post-imputation summaries. |
| **[ADMIXTURE](https://dalexander.github.io/admixture/)** | Model-based ancestry estimation; produces **cross-validation logs** (`cv_log.txt`) and **ancestry proportions** (`.Q` files). |

Exact commands and file names depend on your server pipeline; this repository’s R script expects the named inputs below to already exist in the working directory.

---

## What the R script does

1. **Before imputation** — Reads PLINK `.frq` and `.hwe` (e.g. after MAF / missingness filters), computes He, Ho, MAF, PIC, call rate, summary tables, and four-panel distribution plots.
2. **After imputation** — Same metrics and plots using imputed PLINK outputs.
3. **Choice of *K*** — Plots ADMIXTURE cross-validation error from `cv_log.txt` and **BIC** from `adegenet::find.clusters()` on a `genlight` object built from `imputed.raw`.
4. **ADMIXTURE** — Barplot of *Q* matrix (e.g. `imputed_K4.Q`) with sample order and a 59% assignment rule; writes `admixture_summary.csv`.
5. **UPGMA** — Reads a distance matrix and ID file (e.g. from PLINK `--distance`), builds a fan dendrogram coloured by ADMIXTURE clusters.
6. **DAPC** — Discriminant analysis of principal components (`adegenet`), plot and cluster assignments CSV.
7. **PCA** — PCA on the same `genlight` object (`glPca`), PC1–PC2 plot, scores and variance explained CSVs.

All tabular and figure outputs are written under a **`Results/`** folder (created if missing).

---

## R packages

- `ggplot2`, `gridExtra`, `readr`
- `adegenet`, `ape`

Install in R if needed, e.g. `install.packages(c("ggplot2", "gridExtra", "readr", "ape"))` and `BiocManager::install("adegenet")` or CRAN as appropriate for your setup.

---

## Expected input files (by section)

Adjust paths in `Population_genetics.R` if your names differ (the script sets `setwd()` at the top).

| File(s) | Used for |
|---------|----------|
| `only_chr_clean_maf_09.frq`, `only_chr_clean_maf_09.hwe` | Pre-imputation summaries |
| `imputed.frq`, `imputed.hwe` | Post-imputation summaries |
| `cv_log.txt` | ADMIXTURE cross-validation |
| `imputed.raw` | Genotype matrix for BIC / DAPC / PCA |
| `imputed_K4.Q`, `imputed.fam` | ADMIXTURE barplot and labels (change `K` and filename if you use another *K*) |
| `imputed_dist.dist`, `imputed_dist.dist.id` | UPGMA distance matrix and sample order |

**Note:** Set `N_ind` in the script to your actual sample count (currently `425`). Set `K` for ADMIXTURE plots to match your chosen *K* and `.Q` filename.

---

## How to run

1. Place all required inputs in the project directory (or update `setwd()`).
2. Open R and source the script, or run:  
   `Rscript Population_genetics.R`  
   from that directory.

---

## Outputs (under `Results/`)

Includes CSV summaries, merged frequency/HWE tables, distribution plots (before/after imputation), ADMIXTURE CV plot, BIC plot, ADMIXTURE barplot, UPGMA dendrogram, DAPC plot, PCA plot, PCA scores, variance explained, DAPC cluster assignments, and admixture summary counts.
