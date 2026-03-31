setwd('/Users/mac/Documents/Cursor_Project/GS/Tuber_quality_data/boild_data_folder')

# bWGR: Official CRAN package for Bayesian Whole-Genome Regression
# Published: Xavier et al. (2019) Bioinformatics 36(6):1957-1959
# Authors: Alencar Xavier, William Muir, David Habier, Kyle Kocak, Shizhong Xu, Katy Rainey
# Package provides true BayesCπ implementation (not available in BGLR 1.1.4)
library(bWGR)

# Load data
df_index = read.csv('index_boild_mapping.csv')
df_dosage = read.csv('dosage_boild_converted.csv')

# Define the function to calculate the R² score (coefficient of determination)
r2_score <- function(y_true, y_pred) {
  ss_res <- sum((y_true - y_pred)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  return(1 - (ss_res / ss_tot))
}

# Prepare genotype data
# Set genotype IDs as row names
rownames(df_dosage) <- df_dosage$geno
# Convert to matrix (remove 'geno' column, keep only markers)
X <- as.matrix(df_dosage[, -1])  # Marker matrix: n genotypes × p markers

# Check alignment
if (!all(df_index$geno %in% rownames(X))) {
  stop("Error: Some genotypes in index file are not in dosage file")
}

# Ensure X and df_index are aligned
df_index <- df_index[match(rownames(X), df_index$geno), ]

# Set seed for reproducibility
set.seed(123)

# Define 10-fold cross-validation
n <- nrow(df_index)
k <- 10
folds <- sample(rep(1:k, length.out = n))

# Initialize metrics data frame
matrics <- data.frame(folds = c(1:10), r = NA, r2 = NA, rmse = NA, mae = NA)

# nIter = 50000, burnIn = 10000 for publication-quality results
# nIter = 10000, burnIn = 2000 for quick testing

# Perform 10-fold cross-validation
for (i in 1:k) {
  test_indices <- which(folds == i)
  train_indices <- setdiff(1:n, test_indices)
  
  # Split data into training and test sets
  y_train <- df_index$MGIDI[train_indices]
  X_train <- X[train_indices, ]
  X_test <- X[test_indices, ]
  
  # Fit BayesCπ model using bWGR package
  # BayesCπ: Common variance + π parameter for probability of zero effect (variable selection)
  # bWGR::BayesCpi implements the true BayesCπ model
  # Note: bWGR creates temporary .dat files during MCMC (harmless warnings may appear)
  fm <- suppressWarnings(BayesCpi(y = y_train, 
                                 X = X_train, 
                                 it = 10000,    # nIter: number of MCMC iterations
                                 bi = 2000,     # burnIn: iterations to discard before sampling
                                 df = 5,        # df: degrees of freedom for prior distribution
                                                #      Lower df = more informative prior (tighter)
                                                #      Higher df = less informative prior (wider)
                                                #      Default = 5 (standard in genomic selection)
                                 R2 = 0.5))     # R2: Expected proportion of variance explained by markers
                                                #     Used to set the scale parameter of the prior
                                                #     R2 = 0.5 means we expect markers to explain ~50% of variance
                                                #     This is a common default assumption in genomic selection
                                                #     Can be adjusted based on trait heritability (h²)
                                                #     If h² ≈ 0.3, use R2 = 0.3; if h² ≈ 0.7, use R2 = 0.7
  
  # Extract predictions for test set
  # Predictions = intercept (mu) + X_test %*% coefficients (b)
  predictions <- as.numeric(fm$mu) + X_test %*% fm$b
  predictions <- as.numeric(predictions)  # Convert to vector
  
  # Actual values for test set
  actuals <- df_index$MGIDI[test_indices]
  
  # Calculate prediction accuracy metrics
  r <- cor(actuals, predictions, use = "complete.obs")
  r2 <- r2_score(actuals, predictions)
  rmse <- sqrt(mean((predictions - actuals)^2))
  mae <- mean(abs(predictions - actuals))
  
  # Store metrics
  matrics[i, ] <- c(i, r, r2, rmse, mae)
  
  print(paste("Fold", i, "completed"))
}

# Create Results folder if it doesn't exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Save the metrics to a csv file
write.csv(matrics, "Results/10foldcv_bayescpi.csv", row.names = FALSE)

# Print summary
print("BayesCπ Cross-Validation Results (using bWGR package):")
print(paste("Mean R:", round(mean(matrics$r), 4)))
print(paste("Mean R²:", round(mean(matrics$r2), 4)))
print(paste("Mean RMSE:", round(mean(matrics$rmse), 4)))
print(paste("Mean MAE:", round(mean(matrics$mae), 4)))

# Clean up temporary files created by bWGR during MCMC


