setwd('/Users/mac/Documents/Cursor_Project/GS/Tuber_quality_data/boild_data_folder')

library(BGLR)      # for Bayesian methods

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
  
  # Create phenotype vector with test set masked as NA
  y <- df_index$MGIDI
  y[test_indices] <- NA  # Mask the test set phenotypes
  
  # Fit BayesB model
  # BayesB: Mixture distribution - some markers have zero effect (sparse model)
  fm <- BGLR(y = y, 
             ETA = list(Markers = list(X = X, model = "BayesB")),
             nIter = 10000, 
             burnIn = 2000, 
             verbose = FALSE)
  
  # Extract predictions and actual values for test set
  predictions <- fm$yHat[test_indices]
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
write.csv(matrics, "Results/10foldcv_bayesb.csv", row.names = FALSE)

# Print summary
print("BayesB Cross-Validation Results:")
print(paste("Mean R:", round(mean(matrics$r), 4)))
print(paste("Mean R²:", round(mean(matrics$r2), 4)))
print(paste("Mean RMSE:", round(mean(matrics$rmse), 4)))
print(paste("Mean MAE:", round(mean(matrics$mae), 4)))

