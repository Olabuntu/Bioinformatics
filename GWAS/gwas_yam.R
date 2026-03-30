#
# Yam GWAS workflow (GAPIT). Simulated genotypes + phenotypes (nothing to upload).
# Real study: same steps after PLINK QC on the chip data.
#
# Usage:  cd GWAS && Rscript gwas_yam.R
#

options(warn = 1)

# =============================================================================
# 0) Parameters
# =============================================================================
SEED <- 42
N_LINES <- 300
N_SNPS <- 800
N_CHR <- 10
GAPIT_MODELS <- c("GLM", "FarmCPU")
PCA_TOTAL <- 4
ALPHA <- 0.05
REGION_GAP_BP <- 250000

# QC cutoffs (same idea as PLINK --maf / --geno / --mind / --hwe)
QC_MAF_MIN <- 0.05
QC_SNP_MISS_MAX <- 0.10
QC_IND_MISS_MAX <- 0.10
QC_HWE_P_MIN <- 1e-6

set.seed(SEED)

# =============================================================================
# 1) Packages — install if missing, then load
# =============================================================================
# remotes: install GAPIT from GitHub
# GAPIT:   GWAS + kinship + PCA plots
# utils:   read.csv, write.csv, write.table
# stats:   sd(), pchisq() for HWE

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
library(remotes)

if (!suppressWarnings(requireNamespace("GAPIT", quietly = TRUE))) {
  try(install_github("jiabowang/GAPIT"), silent = TRUE)
}
if (requireNamespace("GAPIT", quietly = TRUE)) {
  library(GAPIT)
} else {
  suppressWarnings(
    source("https://raw.githubusercontent.com/jiabowang/GAPIT/refs/heads/master/gapit_functions.txt", encoding = "UTF-8")
  )
}

library(utils)
library(stats)

# =============================================================================
# 2) Simulate genotype matrix + map + phenotypes
# =============================================================================
line_id <- sprintf("YAM%04d", seq_len(N_LINES))
snp_id <- sprintf("Snp_%05d", seq_len(N_SNPS))

chr <- as.integer(ceiling(seq(from = 0.51, to = N_CHR + 0.49, length.out = N_SNPS)))
pos <- as.integer(sort(sample(1:2e6, N_SNPS)))
map <- data.frame(SNP = snp_id, Chromosome = chr, Position = pos, stringsAsFactors = FALSE)

G <- matrix(sample(0:2, N_LINES * N_SNPS, replace = TRUE), N_LINES, N_SNPS)
colnames(G) <- snp_id
rownames(G) <- line_id

qtn_y <- sort(sample(N_SNPS, 5))
qtn_s <- sort(sample(setdiff(seq_len(N_SNPS), qtn_y), 4))
yield <- 12 + as.vector(G[, qtn_y, drop = FALSE] %*% rnorm(length(qtn_y), 0.35)) + rnorm(N_LINES, 0, 1.1)
stress <- 5 + as.vector(G[, qtn_s, drop = FALSE] %*% rnorm(length(qtn_s), 0.4)) + rnorm(N_LINES, 0, 1.0)
pheno <- data.frame(Taxa = line_id, yield = yield, abiotic_stress = stress, stringsAsFactors = FALSE)

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
write.table(pheno, "data/derived/pheno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(map, "data/derived/map.txt", sep = "\t", quote = FALSE, row.names = FALSE)

meta <- data.frame(
  item = c("lines", "snps", "qtn_index_yield", "qtn_index_stress", "seed"),
  value = c(N_LINES, N_SNPS, paste(qtn_y, collapse = ","), paste(qtn_s, collapse = ","), SEED),
  stringsAsFactors = FALSE
)
write.table(meta, "data/derived/simulation_meta.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# 3) QC (before association) — sample missingness, SNP MAF, SNP missing, HWE
# =============================================================================
# Per-line missing rate across SNPs
miss_ind <- rowMeans(is.na(G) | G < 0 | G > 2)
keep_lines <- miss_ind <= QC_IND_MISS_MAX
n_drop_lines <- sum(!keep_lines)
if (n_drop_lines > 0L) {
  G <- G[keep_lines, , drop = FALSE]
  pheno <- pheno[keep_lines, , drop = FALSE]
  line_id <- pheno$Taxa
}

# Per-SNP missing
miss_snp <- colMeans(is.na(G) | G < 0 | G > 2)
# MAF from dosage (0,1,2)
p_allele <- colMeans(G, na.rm = TRUE) / 2
maf <- pmin(p_allele, 1 - p_allele)

# HWE exact-ish: chi-square test for 3 genotypes at each SNP (biallelic)
hwe_p <- rep(NA_real_, ncol(G))
for (j in seq_len(ncol(G))) {
  g <- G[, j]
  g <- g[!is.na(g)]
  n0 <- sum(g == 0)
  n1 <- sum(g == 1)
  n2 <- sum(g == 2)
  n <- n0 + n1 + n2
  if (n < 5L) {
    next
  }
  p <- (2 * n2 + n1) / (2 * n)
  e0 <- (1 - p)^2 * n
  e1 <- 2 * p * (1 - p) * n
  e2 <- p^2 * n
  if (any(e0 < 1e-6) || any(e1 < 1e-6) || any(e2 < 1e-6)) {
    next
  }
  chisq <- sum((c(n0, n1, n2) - c(e0, e1, e2))^2 / c(e0, e1, e2))
    hwe_p[j] <- pchisq(chisq, df = 1, lower.tail = FALSE)
}

keep_snps <- (maf >= QC_MAF_MIN) &
  (miss_snp <= QC_SNP_MISS_MAX) &
  (is.na(hwe_p) | hwe_p >= QC_HWE_P_MIN)

n_snps_before <- ncol(G)
G <- G[, keep_snps, drop = FALSE]
map <- map[keep_snps, , drop = FALSE]
snp_id <- colnames(G)

geno <- data.frame(Taxa = line_id, G, check.names = FALSE)
write.table(geno, "data/derived/geno.txt", sep = "\t", quote = FALSE, row.names = FALSE)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

qc_snps <- data.frame(
  SNP = snp_id,
  MAF = round(maf[keep_snps], 4),
  missing_rate = miss_snp[keep_snps],
  HWE_p = hwe_p[keep_snps],
  stringsAsFactors = FALSE
)
write.csv(qc_snps, "results/tables/qc_per_snp_after_filters.csv", row.names = FALSE)

qc_report <- c(
  "=== QC report (pre-GWAS) ===",
  paste("Lines at start:", N_LINES, "  dropped (missing):", n_drop_lines, "  kept:", nrow(G)),
  paste("SNPs at start:", n_snps_before, "  dropped (MAF / missing / HWE):", sum(!keep_snps), "  kept:", ncol(G)),
  "",
  "Thresholds used:",
  paste("  MAF >=", QC_MAF_MIN),
  paste("  SNP missing <=", QC_SNP_MISS_MAX),
  paste("  line missing <=", QC_IND_MISS_MAX),
  paste("  HWE p >=", QC_HWE_P_MIN, "(SNPs with stronger departure from HWE removed, like PLINK --hwe)"),
  "",
  "Phenotypes:",
  paste("  yield: mean =", round(mean(pheno$yield), 3), " sd =", round(sd(pheno$yield), 3)),
  paste("  abiotic_stress: mean =", round(mean(pheno$abiotic_stress), 3), " sd =", round(sd(pheno$abiotic_stress), 3)),
  ""
)
writeLines(qc_report, "results/qc_report.txt")
writeLines(qc_report, "results/qc_summary.txt")

# =============================================================================
# 4) GAPIT (association + kinship + PCA inside GAPIT)
# =============================================================================
gapit_dir <- normalizePath("results/gapit", mustWork = FALSE)
dir.create(gapit_dir, recursive = TRUE, showWarnings = FALSE)
old <- getwd()
setwd(gapit_dir)
on.exit(setwd(old), add = TRUE)

traits <- c("yield", "abiotic_stress")
for (tr in traits) {
  message("GAPIT trait: ", tr)
  GAPIT(
    Y = pheno[, c("Taxa", tr), drop = FALSE],
    GD = geno,
    GM = map,
    model = GAPIT_MODELS,
    PCA.total = PCA_TOTAL,
    file.output = TRUE,
    Inter.Plot = FALSE,
    Multiple_analysis = FALSE,
    PCA.3d = FALSE
  )
}

# =============================================================================
# 5) Post-GWAS: Bonferroni, regions, tables
# =============================================================================
region_ids_for_pos <- function(pos, max_gap) {
  ord <- order(pos)
  po <- pos[ord]
  gid <- integer(length(po))
  gid[1] <- 1L
  g <- 1L
  if (length(po) > 1L) {
    for (i in 2:length(po)) {
      if (po[i] - po[i - 1] <= max_gap) {
        gid[i] <- g
      } else {
        g <- g + 1L
        gid[i] <- g
      }
    }
  }
  out <- integer(length(pos))
  out[ord] <- gid
  out
}

gwas_files <- list.files(gapit_dir, pattern = "GWAS_Results", full.names = TRUE)
gwas_files <- gwas_files[grepl("\\.csv$", gwas_files, ignore.case = TRUE)]
gwas_files <- gwas_files[!grepl("StdErr", gwas_files, fixed = TRUE)]
gwas_nyc <- gwas_files[grepl("(NYC)", basename(gwas_files), fixed = TRUE)]
if (length(gwas_nyc) >= 1L) {
  gwas_files <- gwas_nyc
}
if (length(gwas_files) < 1L) {
  stop("No GWAS_Results csv under results/gapit — check GAPIT run.")
}

all_hits <- list()
run_lines <- character()

for (f in gwas_files) {
  tab <- utils::read.csv(f, stringsAsFactors = FALSE)
  pcol <- if ("P.value" %in% names(tab)) "P.value" else if ("P_value" %in% names(tab)) "P_value" else NA_character_
  if (is.na(pcol)) next
  p <- suppressWarnings(as.numeric(tab[[pcol]]))
  m <- sum(is.finite(p))
  thr <- ALPHA / m
  sig <- tab[is.finite(p) & p <= thr, , drop = FALSE]
  sig$p_threshold <- thr
  sig$source_file <- basename(f)
  if (nrow(sig) > 0L) {
    all_hits[[basename(f)]] <- sig
  }

  bn <- basename(f)
  model_guess <- if (grepl("FarmCPU", bn, fixed = TRUE)) "FarmCPU" else "GLM"
  trait_guess <- sub(".*GWAS_Results\\.(?:GLM|FarmCPU)\\.([^.]+)\\(.*", "\\1", bn, perl = TRUE)
  if (trait_guess == bn) trait_guess <- "unknown"

  run_lines <- c(
    run_lines,
    paste(bn, "| tested SNPs:", m, "| Bonferroni alpha:", format(thr, scientific = TRUE),
          "| significant:", nrow(sig))
  )

  if (nrow(sig) < 1L) next

  chr_col <- if ("Chr" %in% names(sig)) "Chr" else "Chromosome"
  pos_col <- if ("Pos" %in% names(sig)) "Pos" else "Position"
  snp_col <- "SNP"
  ps <- suppressWarnings(as.numeric(sig[[pcol]]))
  ch <- suppressWarnings(as.integer(sig[[chr_col]]))
  bp <- suppressWarnings(as.numeric(sig[[pos_col]]))

  leads <- list()
  for (cc in sort(unique(ch))) {
    w <- which(ch == cc)
    gid <- region_ids_for_pos(bp[w], REGION_GAP_BP)
    for (k in unique(gid)) {
      idx <- w[which(gid == k)]
      best <- idx[which.min(ps[idx])]
      leads[[length(leads) + 1L]] <- data.frame(
        trait = trait_guess,
        model = model_guess,
        chromosome = ch[best],
        region_start = min(bp[idx], na.rm = TRUE),
        region_end = max(bp[idx], na.rm = TRUE),
        n_snps = length(idx),
        lead_SNP = as.character(sig[[snp_col]][best]),
        lead_p = ps[best],
        stringsAsFactors = FALSE
      )
    }
  }
  reg <- do.call(rbind, leads)
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", tools::file_path_sans_ext(bn))
  out_reg <- file.path(old, "results/tables", paste0("regions_", safe_name, ".csv"))
  utils::write.csv(reg, out_reg, row.names = FALSE)
}

if (length(all_hits) > 0L) {
  want <- c("SNP", "Chr", "Pos", "P.value", "p_threshold", "source_file")
  trimmed <- lapply(all_hits, function(z) {
    cn <- names(z)
    if (!"P.value" %in% cn && "P_value" %in% cn) {
      names(z)[names(z) == "P_value"] <- "P.value"
    }
    z[, intersect(want, names(z)), drop = FALSE]
  })
  combined <- do.call(rbind, trimmed)
  utils::write.csv(combined, file.path(old, "results/tables/significant_snps_all_models.csv"), row.names = FALSE)
}

pca_blurb <- character()
pca_file <- file.path(gapit_dir, "GAPIT.Genotype.PCA_eigenvalues.csv")
if (file.exists(pca_file)) {
  pe <- tryCatch(utils::read.csv(pca_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(pe) && nrow(pe) > 0) {
    pca_blurb <- c(
      "PCA (first axes, from GAPIT)",
      "-------------------------------------------",
      paste(utils::capture.output(print(pe[seq_len(min(5L, nrow(pe))), , drop = FALSE])), collapse = "\n"),
      ""
    )
  }
}

notes <- c(
  "What this repository shows",
  "-------------------------",
  "Simulated genotypes/phenotypes so nothing confidential is uploaded.",
  "Same file layout and QC rules as the real panel after PLINK filtering.",
  "",
  "Outputs",
  "-------",
  "  results/qc_report.txt        QC before GWAS",
  "  results/gapit/               GAPIT plots and GWAS csv",
  "  results/tables/              filtered SNP list, significant hits, regions",
  "",
  paste("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  ""
)
writeLines(c(notes, pca_blurb, "Per-file significance:", run_lines), file.path(old, "results/STUDY_NOTES.txt"))

message("Done. See results/qc_report.txt and results/gapit/")
