setwd('/Users/mac/Documents/Cursor_Project/GS/Tuber_quality_data/boild_data_folder')

library(rrBLUP)   # for A.mat()
library(BGLR)      # for BGLR() and RKHS model



df_index = read.csv('index_boild_mapping.csv')
df_dosage = read.csv('dosage_boild_converted.csv')



# define the function to calculate the R² score
r2_score <- function(y_true, y_pred) {
  ss_res <- sum((y_true - y_pred)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  return(1 - (ss_res / ss_tot))
}




# performing additive genomic relationship matrix (G) estimation for GS analysis

rownames(df_dosage) <- df_dosage$geno
df_dosage <- as.matrix(df_dosage[, -c(1)])  # Exclude the first column (geno)   # convert the dataframe to a matrix

if (!all(df_index$geno %in% rownames(df_dosage))) {
  stop("Error: Some genotypes in index file are not in dosage file")
}



# Perform the additive genomic relationship matrix (G) estimation
A <- A.mat(df_dosage)

# ============================================================================
# G×E MODELING REQUIREMENTS CHECK
# ============================================================================
# For G×E modeling, you need phenotypic data with:
# - Multiple environments (ENV column)
# - Multiple observations per genotype (same geno in different ENV)
# - Year information (optional but recommended)
#
# Your current data (df_index) has only: geno, MGIDI
# This is AGGREGATED data (one value per genotype) - NOT suitable for G×E
#
# Required data structure:
#   geno        | ENV  | Year | MGIDI
#   TDr1100873  | ENV1 | 2021 | 2.3
#   TDr1100873  | ENV2 | 2021 | 2.5
#   TDr1100873  | ENV1 | 2022 | 2.4
#   ...
# ============================================================================

# Check if data has environment structure
if (!("ENV" %in% colnames(df_index)) && !("Year" %in% colnames(df_index))) {
  stop("ERROR: G×E modeling requires ENV and/or Year columns in phenotypic data.\n",
       "Your current data has only geno and MGIDI (aggregated values).\n",
       "G×E cannot be modeled with this data structure.\n",
       "Please use: additive.R, additive_dominance.R, or Additive_epistatic.R instead.")
}

# If data has environment structure, proceed with G×E setup
if ("ENV" %in% colnames(df_index) || "Year" %in% colnames(df_index)) {
  # Create environment factor
  if ("ENV" %in% colnames(df_index) && "Year" %in% colnames(df_index)) {
    df_index$ENV_factor <- as.factor(paste(df_index$ENV, df_index$Year, sep = "_"))
  } else if ("ENV" %in% colnames(df_index)) {
    df_index$ENV_factor <- as.factor(df_index$ENV)
  } else {
    df_index$ENV_factor <- as.factor(df_index$Year)
  }
  
  unique_env <- unique(df_index$ENV_factor)
  n_env <- length(unique_env)
  
  # Create environment-specific genomic relationship matrices for G×E
  # Each environment gets its own kernel matrix (using same A matrix)
  GXE <- vector("list", n_env)
  names(GXE) <- unique_env
  
  for (e in unique_env) {
    # For G×E modeling, we use the same genomic relationship matrix
    # but BGLR will estimate separate variance components for each environment
    GXE[[e]] <- A  
  }
  
  # Define the ETA list for BGLR with additive + G×E effects
  ETA <- list(Add = list(K = A, model = "RKHS"))
  
  # Add environment-specific G×E terms
  for (e in unique_env) {
    ETA[[paste0("GxE_", e)]] <- list(K = GXE[[e]], model = "RKHS")
  }
} else {
  # Fallback: Use additive only if G×E not possible
  warning("G×E modeling not possible. Using additive model only.")
  ETA <- list(Add = list(K = A, model = "RKHS"))
}

set.seed(123)

# Define 10 folds
n <- dim(df_index)[1]
k <- 10
folds <- sample(rep(1:k, length.out = n))

matrics <- data.frame( folds = c(1:10) , r = NA, r2 = NA, rmse = NA, mae = NA)

# Perform 10-fold cross-validation
for (i in 1:k) {
  test_indices <- which(folds == i)
  
  # Create phenotype vector with test set masked as NA
  y_train <- df_index$MGIDI
  y_train[test_indices] <- NA
  
  # nIter = 50000, burnIn = 10000 for publication-quality results
  # nIter = 10000, burnIn = 2000 for quick testing
  
  # Fit BGLR model with additive + G×E effects
  fm <- BGLR(y = y_train, ETA = ETA, nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  # Extract predictions and actual values for test set
  preds <- fm$yHat[test_indices] 
  obs <- df_index$MGIDI[test_indices]
  
  # Calculate prediction accuracy metrics
  r <- cor(obs, preds, use = "complete.obs")
  r2 <- r2_score(obs, preds)  # Correct order: true first, then predicted
  # calculate the RMSE
  rmse <- sqrt(mean((preds - obs)^2))
  # calculate the MAE
  mae <- mean(abs(preds - obs))   
  
  # Store metrics
  matrics[i, ] <- c(i, r, r2, rmse, mae)
  
  print(paste("Fold", i, "completed"))
}

# Create Results folder if it doesn't exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Save cross-validation results
write.csv(matrics, "Results/10foldcv_additive_GXE_gblup.csv", row.names = FALSE)