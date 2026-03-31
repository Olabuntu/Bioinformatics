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

# Calculate epistatic relationship matrix using Hadamard product approach
# This captures pairwise epistatic interactions between loci
# E = A ⊙ A (Hadamard product) captures epistatic relationships
E <- A * A  # Element-wise multiplication for epistatic relationship matrix

# Normalize epistatic matrix to have mean diagonal of 1
# This ensures the matrix is on a similar scale to the additive matrix
diag_mean <- sum(diag(E)) / nrow(E)
if (diag_mean > 0) {
  E <- E / diag_mean
} else {
  stop("Error: Epistatic matrix has zero diagonal mean")
}



set.seed(123)

# Define 10 folds
n <- dim(df_index)[1]
k <- 10
folds <- sample(rep(1:k, length.out = n))


matrics <- data.frame( folds = c(1:10) , r = NA, r2 = NA, rmse = NA, mae = NA)

# Define the ETA list with both additive and epistatic effects
# Using RKHS (Reproducing Kernel Hilbert Space) model for both
ETA <- list(Add = list(K = A, model = 'RKHS'),
            Epi = list(K = E, model = 'RKHS'))


# nIter = 50000, burnIn = 10000 for publication-quality results
# nIter = 10000, burnIn = 2000 for quick testing

# Perform 10-fold cross-validation
for (i in 1:k) {
  
  test_indices <- which(folds == i)
  
  # Create phenotype vector with test set masked as NA
  y <- df_index$MGIDI
  y[test_indices] <- NA  # Mask the test set phenotypes
  
  # Fit BGLR model with additive and epistatic effects
  fm <- BGLR(y = y, ETA = ETA, nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  # Extract predictions and actual values for test set
  predictions <- fm$yHat[test_indices]
  actuals <- df_index$MGIDI[test_indices]
  
  # Calculate prediction accuracy metrics
  r <- cor(actuals, predictions, use = "complete.obs")
  r2 <- r2_score(actuals, predictions)
  # calculate the RMSE
  rmse <- sqrt(mean((predictions - actuals)^2))
  # calculate the MAE
  mae <- mean(abs(predictions - actuals))
  
  # Store metrics
  matrics[i, ] <- c(i, r, r2, rmse, mae)
  
  print(paste("Fold", i, "completed"))
}

# Create Results folder if it doesn't exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# save the metrics to a csv file
write.csv(matrics, "Results/10foldcv_additive_epistatic_gblup.csv", row.names = FALSE) # save the metrics to a csv file
