setwd('/Users/mac/Documents/Cursor_Project/GS/Tuber_quality_data/boild_data_folder')

library(BGLR)      # for BGLR() and RKHS model
library(AGHmatrix) # for Gmatrix() with method = "Vitezica"



df_index = read.csv('index_boild_mapping.csv')
df_dosage = read.csv('dosage_boild_converted.csv')



# define the function to calculate the R² score
r2_score <- function(y_true, y_pred) {
  ss_res <- sum((y_true - y_pred)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  return(1 - (ss_res / ss_tot))
}




# Performing dominance-only genomic relationship matrix (D) estimation for GS analysis

rownames(df_dosage) <- df_dosage$geno
df_dosage <- as.matrix(df_dosage[, -c(1)])  # Exclude the first column (geno)   # convert the dataframe to a matrix

if (!all(df_index$geno %in% rownames(df_dosage))) {
  stop("Error: Some genotypes in index file are not in dosage file")
}




# Perform the dominance genomic relationship matrix (D) estimation
# Using Vitezica method for calculating dominance relationships
D <- Gmatrix(df_dosage, method = "Vitezica", ploidy = 2)

# Define the ETA list with dominance effects only
ETA <- list(Dom = list(K = D, model = "RKHS"))

# set the seed for reproducibility
set.seed(123)

# define the number of folds
n <- dim(df_index)[1]
k <- 10 # 10 folds
folds <- sample(rep(1:k, length.out = n))


# initialize the metrics data frame
matrics <- data.frame( folds = c(1:10) , r = NA, r2 = NA, rmse = NA, mae = NA)


# Perform 10-fold cross-validation
for (i in 1:k) {
  
  test_indices <- which(folds == i)
  
  # Create phenotype vector with test set masked as NA
  y <- df_index$MGIDI
  y[test_indices] <- NA  # Mask the test set phenotypes
  
  # nIter = 50000, burnIn = 10000 for publication-quality results
  # nIter = 10000, burnIn = 2000 for quick testing
  
  # Fit BGLR model with dominance effects only
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

# Save the metrics to a csv file
write.csv(matrics, "Results/10foldcv_dominance_gblup.csv", row.names = FALSE)
