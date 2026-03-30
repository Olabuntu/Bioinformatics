setwd('/Users/mac/Documents/IITA/Marta_mozanbique/N')

#create a Results folder 
if (!dir.exists("Results")) {
  dir.create("Results")
}


######SUMMARY BEFORE IMPUTATION AFTER FILTERING MAF 0.05 AND MAX MISSINGNESS 0.9 ########
# Summary statistics from PLINK .frq and .hwe (no .raw needed)
library(ggplot2)
library(gridExtra) # for arranging multiple plots
library(readr) # for read_table()
library(adegenet) # for genlight object
library(ape) # for hierarchical clustering


frq <- read_table(
  "only_chr_clean_maf_09.frq",
  col_names = c("CHROM", "POS", "N_ALLELES", "N_CHR", 
                "ALLELE1_FREQ", "ALLELE2_FREQ"),
  skip = 1
)


# 2) Read PLINK Hardy-Weinberg output for observed heterozygotes
hwe <- read_table("only_chr_clean_maf_09.hwe")

p = as.numeric(sub("^[ATGC]:", "", frq$ALLELE1_FREQ))
q = as.numeric(sub("^[ATGC]:", "", frq$ALLELE2_FREQ))

MAF <- pmin(p, q)

He <- 2 * p * q

PIC <- 2 * p * q * (1 - p * q)

#call rate per SNP
N_ind <- 425  # set manually based on your sample size
call_rate <-  ((frq$N_CHR / 2) / N_ind) * 100

# Split the OBS column (name may vary: "OBS(HOM1/HET/HOM2)" or similar)
obs_col <- names(hwe)[grep("OBS.*HOM.*HET", names(hwe))][1]
obs_split <- strsplit(hwe[[obs_col]], "/", fixed = TRUE)

# Convert to numeric and compute Ho
obs_mat <- do.call(rbind, lapply(obs_split, function(x) as.numeric(x)))
HOM1 <- obs_mat[, 1]
HET  <- obs_mat[, 2]
HOM2 <- obs_mat[, 3]
N_obs <- HOM1 + HET + HOM2

Ho <- HET / N_obs   # observed proportion heterozygous



# Merge frq and hwe 
frq_df <- data.frame(
  CHROM     = frq$CHROM,
  POS       = frq$POS,
  N_CHR     = frq$N_CHR,
  MAF       = MAF,
  He        = He,
  PIC       = PIC,
  call_rate = call_rate,
  stringsAsFactors = FALSE
)
hwe_df <- data.frame(
  CHROM = hwe$CHR,
  POS   = hwe$POS,
  Ho    = Ho,
  P_HWE = hwe$P_HWE,
  stringsAsFactors = FALSE
)
merged <- merge(frq_df, hwe_df, by = c("CHROM", "POS"), all.x = TRUE)



#  Summary statistics 
summary_stats <- data.frame(
  Statistic = c("He", "Ho", "MAF", "PIC", "Call Rate (%)"),
  Mean = c(mean(merged$He, na.rm = TRUE),
           mean(merged$Ho, na.rm = TRUE),
           mean(merged$MAF, na.rm = TRUE),
           mean(merged$PIC, na.rm = TRUE),
           mean(merged$call_rate, na.rm = TRUE)),
  Min = c(min(merged$He, na.rm = TRUE),
          min(merged$Ho, na.rm = TRUE),
          min(merged$MAF, na.rm = TRUE),
          min(merged$PIC, na.rm = TRUE),
          min(merged$call_rate, na.rm = TRUE)),
  Max = c(max(merged$He, na.rm = TRUE),
          max(merged$Ho, na.rm = TRUE),
          max(merged$MAF, na.rm = TRUE),
          max(merged$PIC, na.rm = TRUE),
          max(merged$call_rate, na.rm = TRUE))
)


# Fig. 4-style: four panels (He, Ho, MAF, PIC) ---
# A: He – green, no border
p_he <- ggplot(merged, aes(x = He)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "darkgreen", alpha = 0.8, color = "black", linewidth = 0.5) +
  labs(title = "Distribution of expected heterozygosity (He)", x = "Expected heterozygosity", y = "Density")

# B: Ho – red, red border
p_ho <- ggplot(merged, aes(x = Ho)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "red", alpha = 0.8, color = "red", linewidth = 0.5) +
  labs(title = "Distribution of observed heterozygosity (Ho)", x = "Observed heterozygosity", y = "Density")

# C & D: MAF, PIC – teal with black border
p_maf <- ggplot(merged, aes(x = MAF)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "cadetblue3", alpha = 0.8, color = "black", linewidth = 0.5) +
  labs(title = "Distribution of minor allele frequency", x = "Minor allele frequency", y = "Density of SNPs")

p_pic <- ggplot(merged, aes(x = PIC)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "cadetblue3", alpha = 0.8, color = "black", linewidth = 0.5) +
  labs(title = "Distribution of polymorphic information content", x = "PIC", y = "Density of SNPs")

grid <- grid.arrange(p_he, p_ho, p_maf, p_pic, ncol = 2)


#save summary statistics, and the merged data and the plots
write.csv(summary_stats, "Results/summary_statistics_before_imputation.csv", row.names = FALSE)
write.csv(merged, "Results/merged_frq_hwe_before_imputation.csv", row.names = FALSE)
ggsave("Results/distribution_plots_before_imputation.png", grid, width = 12, height = 10, dpi = 300)






#### SUMMARY AFTER IMPUTATION WITH BEAGLE ########

frq_imp <- read_table(
  "imputed.frq",
  col_names = c("CHROM", "POS", "N_ALLELES", "N_CHR", 
                "ALLELE1_FREQ", "ALLELE2_FREQ"),
  skip = 1
)

hwe_imp <- read_table("imputed.hwe")

p_imp = as.numeric(sub("^[ATGC]:", "", frq_imp$ALLELE1_FREQ))
q_imp = as.numeric(sub("^[ATGC]:", "", frq_imp$ALLELE2_FREQ))

MAF_imp <- pmin(p_imp, q_imp)
He_imp <- 2 * p_imp * q_imp
PIC_imp <- 2 * p_imp * q_imp * (1 - p_imp * q_imp)

#call rate per SNP
N_ind <- 425  # set manually based on your sample size
call_rate_imp <-  ((frq_imp$N_CHR / 2) / N_ind) * 100

# Split the OBS column (name may vary: "OBS(HOM1/HET/HOM2)" or similar)
obs_col_imp <- names(hwe_imp)[grep("OBS.*HOM.*HET", names(hwe_imp))][1]
obs_split_imp <- strsplit(hwe_imp[[obs_col_imp]], "/", fixed = TRUE)

# Convert to numeric and compute Ho
obs_mat_imp <- do.call(rbind, lapply(obs_split_imp, function(x) as.numeric(x)))
HOM1_imp <- obs_mat_imp[, 1]
HET_imp  <- obs_mat_imp[, 2]
HOM2_imp <- obs_mat_imp[, 3]       
N_obs_imp <- HOM1_imp + HET_imp + HOM2_imp
Ho_imp <- HET_imp / N_obs_imp   # observed proportion heterozygous

# Merge frq and hwe
frq_df_imp <- data.frame(
  CHROM     = frq_imp$CHROM,
  POS       = frq_imp$POS,
  N_CHR     = frq_imp$N_CHR,
  MAF       = MAF_imp,
  He        = He_imp,
  PIC       = PIC_imp,
  call_rate = call_rate_imp,
  stringsAsFactors = FALSE
)

hwe_df_imp <- data.frame(
  CHROM = hwe_imp$CHR,
  POS   = hwe_imp$POS,
  Ho    = Ho_imp,
  P_HWE = hwe_imp$P_HWE,
  stringsAsFactors = FALSE
)


merged_imp <- merge(frq_df_imp, hwe_df_imp, by = c("CHROM", "POS"), all.x = TRUE)

# Summary statistics after imputation
summary_stats_imp <- data.frame( 
  Statistic = c("He", "Ho", "MAF", "PIC", "Call Rate (%)"),
  Mean = c(mean(merged_imp$He, na.rm = TRUE),
           mean(merged_imp$Ho, na.rm = TRUE),
           mean(merged_imp$MAF, na.rm = TRUE),
           mean(merged_imp$PIC, na.rm = TRUE),
           mean(merged_imp$call_rate, na.rm = TRUE)),
  Min = c(min(merged_imp$He, na.rm = TRUE),
          min(merged_imp$Ho, na.rm = TRUE),
          min(merged_imp$MAF, na.rm = TRUE),
          min(merged_imp$PIC, na.rm = TRUE),
          min(merged_imp$call_rate, na.rm = TRUE)),
  Max = c(max(merged_imp$He, na.rm = TRUE),
          max(merged_imp$Ho, na.rm = TRUE),
          max(merged_imp$MAF, na.rm = TRUE),
          max(merged_imp$PIC, na.rm = TRUE),
          max(merged_imp$call_rate, na.rm = TRUE))
)

#fig 4-style plots after imputation
p_he_imp <- ggplot(merged_imp, aes(x = He)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "darkgreen", alpha = 0.8, color = "black", linewidth = 0.5) +
  labs(title = "Distribution of expected heterozygosity (He) ", x = "Expected heterozygosity", y = "Density")

p_ho_imp <- ggplot(merged_imp, aes(x = Ho)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "red", alpha = 0.8, color = "red", linewidth = 0.5) +
  labs(title = "Distribution of observed heterozygosity (Ho)", x = "Observed heterozygosity", y = "Density")

p_maf_imp <- ggplot(merged_imp, aes(x = MAF)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "cadetblue3", alpha = 0.8, color = "black", linewidth = 0.5) +
  labs(title = "Distribution of minor allele frequency", x = "Minor allele frequency", y = "Density of SNPs")

p_pic_imp <- ggplot(merged_imp, aes(x = PIC)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "cadetblue3", alpha = 0.8, color = "black", linewidth = 0.5) +
  labs(title = "Distribution of polymorphic information content (PIC) ", x = "PIC", y = "Density of SNPs")
grid_imp <- grid.arrange(p_he_imp, p_ho_imp, p_maf_imp, p_pic_imp, ncol = 2)



# Save summary statistics, merged data, and plots after imputation
write.csv(summary_stats_imp, "Results/summary_statistics_after_imputation.csv", row.names = FALSE)
write.csv(merged_imp, "Results/merged_frq_hwe_after_imputation.csv", row.names = FALSE)
ggsave("Results/distribution_plots_after_imputation.png", grid_imp, width = 12, height = 10, dpi = 300)






######## we need to determine K  from cv for Admixture analysis and from BIC for other analysis


# Your CV errors (K = 1 to 50)

lines <- readLines("cv_log.txt")

# Extract K and CV error with regex
K_adm <- as.numeric(sub(".*K=([0-9]+).*", "\\1", lines))
CV_err <- as.numeric(sub(".*:\\s*([0-9.]+)\\s*$", "\\1", lines))

# Put into a data frame and sort by K
cv_df <- data.frame(K = K_adm, CV_error = CV_err)
cv_df <- cv_df[order(cv_df$K), ]

cv_df

# Data frame for ggplot



# Line + points so the elbow is easy to see
k_admix_plot <- ggplot(cv_df, aes(x = K, y = CV_error)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1.2) +
  scale_x_continuous(breaks = seq(0, 50, by = 5)) +  # or breaks = 1:50 for every K
  labs(
    x = "Number of clusters (K)",
    y = "Cross-validation error",
    title = "CV error vs Clusters (K)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("Results/Admixture_CV_error_plot.png", k_admix_plot, width = 12, height = 8, dpi = 300)



################################## BIC 


# Read: skip first 6 columns (FID, IID, PAT, MAT, SEX, PHENOTYPE)
dat <- read.table('imputed.raw', header = TRUE, stringsAsFactors = FALSE)


# Genotype matrix: individuals x SNPs (only 0/1/2 columns)
geno <- as.matrix(dat[, -(1:6)])
rownames(geno) <- dat$IID   

# Replace NA if coded as -9 or something else in PLINK
any( geno == -9 )  # check if -9 is present)
geno[geno == -9] <- NA # which we don't have after imputation, but just in case


# genlight expects rows = individuals, cols = SNPs, values 0/1/2
x <- new("genlight", geno, ploidy = 2)
#  add SNP names
locNames(x) <- colnames(geno)
indNames(x) <- rownames(geno)


set.seed(42)
clust <- find.clusters(x, max.n.clust = 50, n.pca = NULL, choose.n.clust = FALSE)
# i was asked to pick pc to retain 

bic_by_k <- clust$Kstat
K_vals <- as.numeric(sub("K=", "", names(bic_by_k)))


Bic_df <- data.frame(K = K_vals, BIC = as.numeric(bic_by_k))
bic_plt <- ggplot(Bic_df, aes(x = K, y = BIC)) +
  geom_line(linewidth = 0.8, color = "darkred") +
  geom_point(size = 3, color = "darkred") +
  scale_x_continuous(breaks = K_vals) +
  labs(x = "Number of clusters (K)", y = "BIC",
       title = "BIC vs number of clusters ") +
  theme_minimal()
ggsave("Results/BIC_plot.png", bic_plt, width = 12, height = 6, dpi = 300)


############################################# Admixture plot 
K <- 4  # set K based on CV error plot (elbow method) plot

#  Read Q and sample IDs
Q <- as.matrix(read.table("imputed_K4.Q", header = FALSE))
fam <- read.table("imputed.fam", header = FALSE, stringsAsFactors = FALSE)
ids   <- fam[, 2]
rownames(Q) <- ids
colnames(Q) <- paste0("Cluster", 1:K)

# 59% rule: assigned vs admixed (as in paper)
max_anc    <- apply(Q, 1, max)
assigned   <- max_anc >= 0.59
cluster    <- max.col(Q)
cluster[!assigned] <- NA   # admixed = no single cluster (or use 5 for "Admixed")




ord_by_clust <- list()
for (k in 1:K) {
  idx <- which(assigned & cluster == k)
  if (length(idx) > 0) {
    # Descending order of ancestry in cluster k
    idx <- idx[order(Q[idx, k], decreasing = TRUE)]
    ord_by_clust[[k]] <- idx
  }
}
ord_clust <- unlist(ord_by_clust)
ord_admixed <- which(!assigned)
ord <- c(ord_clust, ord_admixed)

Q_plot <- Q[ord, ]
rownames(Q_plot) <- ids[ord]   # keep sample names for x-axis
cols <- c("darkblue", "green3", "cyan3", "purple")  # one color per cluster


## Save the admixture plot as PNG
png(paste0( "Results/Admixture_plot_K", K, ".png"), width = 18, height = 9, units = "in", res = 300)
# Barplot (leave right margin for legend)
barplot(t(Q_plot), col = cols, border = NA, space = 0,
        xlab = "", ylab = "Ancestry proportion",
        xaxt = "n")   # no default x-axis yet

# Sample names on x-axis: one per bar, vertical, small
n_bars <- nrow(Q_plot)
axis(1, at = (1:n_bars), labels = rownames(Q_plot), 
     tick = FALSE, las = 2, cex.axis = 0.4, line = -0.9)
legend("topright", legend = colnames(Q_plot), fill = cols, bty = "o", 
       inset = -0.002, cex = 0.5)

dev.off()

# save some admixture results for repoeporting
admixture_lower_than_59 <- sum(max_anc < 0.59)
pure_ind_higher_eqal_59 <- sum(max_anc >= 0.59)
admixture_higest_membership_sub_population <- colnames(Q)[max.col(Q)]
max_anc <- apply(Q, 1, max)
higest_value_whole_population <- max(max_anc)
cluster_1 <- sum(admixture_higest_membership_sub_population == "Cluster1")
cluster_2 <- sum(admixture_higest_membership_sub_population == "Cluster2")
cluster_3 <- sum(admixture_higest_membership_sub_population == "Cluster3")
cluster_4 <- sum(admixture_higest_membership_sub_population == "Cluster4")

admixture_summary <- data.frame(
  Cluster = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed", "Pure (>=59%)", "Highest membership in whole population"),
  Count = c(cluster_1, cluster_2, cluster_3, cluster_4, admixture_lower_than_59, pure_ind_higher_eqal_59, higest_value_whole_population)
)

write.csv(admixture_summary, "Results/admixture_summary.csv", row.names = FALSE)



################################# UPGMA dendrogram

# Paths (adjust if needed)
dist_file <- "imputed_dist.dist"
id_file   <- "imputed_dist.dist.id"

# Sample IDs (order = order of rows/columns in the distance matrix)
ids <- read.table(id_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
n   <- length(ids)

# Square matrix: one row per sample, spaces/tabs between values
dist_mat <- as.matrix(read.table(dist_file, header = FALSE))
rownames(dist_mat) <- ids
colnames(dist_mat) <- ids

# Convert to "dist" object for clustering
dist_obj <- as.dist(dist_mat)
# UPGMA hierarchical clustering
upgma <- hclust(dist_obj, method = "average")



Q <- as.matrix(read.table("imputed_K4.Q", header = FALSE))
max_anc <- apply(Q, 1, max)
cluster <- max.col(Q)
cluster[max_anc < 0.59] <- 5L   # 5 = admixed
names(cluster) <- ids


# Convert hclust to phylo and match tip order
tree <- as.phylo(upgma)
tip_colors <- c("darkblue", "green3", "cyan3", "purple", "red")[cluster[tree$tip.label]] # cluster[tree$tip.label] gives the cluster assignment for each tip in the order they appear in the tree, and we use that to index into the color vector. The colors correspond to Cluster 1, Cluster 2, Cluster 3, Cluster 4, and Admixed (5).

# Plot dendrogram
png("Results/UPGMA_dendrogram.png", width = 12, height = 8, units = "in", res = 300)

plot(tree, type = "fan", tip.color = tip_colors, cex = 0.25, label.offset = 0.5)

legend(title = 'Class',"bottomleft", legend = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Admixed"),
       fill = c("darkblue", "green3", "cyan3", "purple", "red"), bty = "n")
dev.off()





###### DAPC (Discriminant Analysis of Principal Components) anlysis

# DAPC requires a genlight object and a grouping factor (clusters)
# We already have the genlight object "x" and the cluster assignments in "cluster" from BIC analysis

set.seed(42)
grp  <- find.clusters(x, n.pca = 150, n.clust = 4, choose.n.clust = FALSE) # picked 150 pcss cause i picked 150 for BIC plot, but you can adjust based on that plot, k = bic k
dapc_res <- dapc(x, grp$grp, n.pca = 150, n.da = 3 ) # n.da = number of discriminant functions to retain (max = K-1)

# plot

cols_dapc <- c("darkblue", "green3", "cyan3", "purple")
png("Results/DAPC_plot.png", width = 12, height = 8, units = "in", res = 300)
scatter(dapc_res, col = cols_dapc, pch = 19, cex = 1.2, legend = TRUE,
        posi.leg = "topright", clab = 0, cstar = 0, scree.da = TRUE,
        bg = "white")

dev.off()

# save cluster assignments for DAPC
dapc_clusters <- data.frame(SampleID = names(grp$grp), DAPC_Cluster = grp$grp)
write.csv(dapc_clusters, "Results/DAPC_cluster_assignments.csv", row.names = FALSE )







########################### PCA
# USING x genlight object from before 

pca_res <- glPca(x, nf = 10, parallel = TRUE, n.cores = 4) # nf = number of principal components to retain

# using the same cluster assignments from BIC for coloring grp$grp

pca_scores <- as.data.frame(pca_res$scores)
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))
pca_scores$IID <- rownames(pca_scores)
pca_scores$Cluster <- grp$grp[match(pca_scores$IID, names(grp$grp))]

# variance explained by each PC
var_explained <- pca_res$eig / sum(pca_res$eig) * 100
var_df <- data.frame(PC = paste0("PC", 1:length(var_explained)), Variance_Explained = var_explained)

# Plot PC1 vs PC2 colored by cluster
cluster_colors <- c("1" = "darkblue", "2" = "green3", "3" = "cyan3", "4" = "purple")
cluster_shapes <- c("1" = 16, "2" = 17, "3" = 15, "4" = 4)
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) +
  geom_point(size = 2.5, alpha = 0.85, stroke = 0.5) +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  scale_shape_manual(values = cluster_shapes, name = "Cluster") +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
    title = "Individual PCA"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


ggsave("Results/PCA_plot.png", p, width = 10, height = 7, dpi = 300)

# Save PCA scores and variance explained
write.csv(pca_scores, "Results/PCA_scores.csv", row.names = FALSE)
write.csv(var_df, "Results/PCA_variance_explained.csv", row.names = FALSE)